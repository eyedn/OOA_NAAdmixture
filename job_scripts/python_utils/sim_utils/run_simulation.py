###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           run_simulation.py
###############################################################################


import pickle
import stdpopsim
import msprime
import tskit
from python_utils import log_msg
from .build_demography import build_demography
from .build_metadata import build_metadata
from python_utils.ancestry_utils.build_global_ancestry_table import (
    build_global_ancestry_table
)
from python_utils.ancestry_utils.get_local_ancestry_table import (
    get_local_ancestry_table
)
from python_utils.ancestry_utils.write_global_ancestry import (
    write_global_ancestry
)
from python_utils.ancestry_utils.write_local_ancestry import (
    write_local_ancestry
)


# simulate the OOA_NAAdmixture demography and generate necessary data files
def run_simulation(args):
    # define species and chromosome
    species = stdpopsim.get_species("HomSap")
    contig = species.get_contig(
        f"chr{args.chromosome}", genetic_map=args.genetic_map
        )
    
    # define sample set for msprime simulator later
    sample_sets = [
        msprime.SampleSet(
            args.sample_size, population="AFR", time=0, ploidy=2
            ),
        msprime.SampleSet(
            args.sample_size, population="EUR", time=0, ploidy=2
            ),
        msprime.SampleSet(
            args.sample_size, population="ADX", time=0, ploidy=2
            ),
    ]

    # define a vector of individual sample names for VCF output
    sample_names = []
    for pop in args.pops:
        for sample_idx in range(args.sample_size):
            sample_names.append(f"{pop}_{sample_idx + 1}")

    # generate demography object
    def parse_float_list(raw_values):
        return [
            float(value) for value in raw_values.split(",") if value != ""
            ]
    demography, model_metadata = build_demography(
        generation_time=args.generation_time,
        mutation_rate=args.mutation_rate,
        t_af_years=args.t_af_years,
        t_ooa_years=args.t_ooa_years,
        t_eu0_years=args.t_eu0_years,
        t_eg_years=args.t_eg_years,
        r_eu0=args.r_eu0,
        r_eu=args.r_eu,
        r_af=args.r_af,
        n_a=args.n_a,
        n_af1=args.n_af1,
        n_b=args.n_b,
        n_eu0=args.n_eu0,
        m_af_b=args.m_af_b,
        m_af_eu=args.m_af_eu,
        admixture_time=args.admixture_time,
        admix_generation_count=args.admix_generation_count,
        admix_mixing_generation_count=args.admix_mixing_generation_count,
        admix_ne_by_generation=parse_float_list(
            args.admix_ne_by_generation
        ),
        admix_afr_props_by_generation=parse_float_list(
            args.admix_afr_props_by_generation
        ),
        admix_eur_props_by_generation=parse_float_list(
            args.admix_eur_props_by_generation
        ),
        admix_prioradmix_props_by_generation=parse_float_list(
            args.admix_prioradmix_props_by_generation
        ),
        admix_modern_growth_rate=args.admix_modern_growth_rate,
        census_time_offset=args.census_time_offset,
    )

    # simulate the above demography with msprime and output the tree sequence
    log_msg("simulating tree sequence")
    ts = msprime.sim_ancestry(
        samples=sample_sets,
        demography=demography,
        recombination_rate=contig.recombination_map,
        sequence_length=contig.length,
        ploidy=2,
        random_seed=args.seed,
        model=args.msprime_model,
    )
    ts = msprime.sim_mutations(
        ts,
        rate=args.mutation_rate,
        random_seed=args.seed + 1,
    )
    if not isinstance(ts, tskit.TreeSequence):
        raise TypeError("Expected a tskit.TreeSequence result")
    ts_path = f"{args.tree_prefix}.ts"
    ts.dump(ts_path)

    # pickle demography and model metadata objects
    log_msg("writing demography and metadata pickles")
    with open(f"{args.pickle_prefix}.demography.pkl", "wb") as out_file:
        pickle.dump(demography, out_file)
    with open(f"{args.pickle_prefix}.model_metadata.pkl", "wb") as out_file:
        pickle.dump(model_metadata, out_file)

    # output vcf resulting from tree sequence
    log_msg("writing all-sample VCF")
    with open(args.vcf_path, "w", encoding="utf-8") as out_file:
        ts.write_vcf(
            out_file,
            contig_id=f"chr{args.chromosome}",
            individual_names=sample_names,
            position_transform=lambda positions: [pos + 1 for pos in positions],
        )

    # output metadata for downstream "admixture --supervised"
    log_msg("writing sample metadata")
    sample_metadata_rows = build_metadata(
        args.pops,
        args.sample_size,
    )
    with open(args.sample_metadata_path, "w", encoding="utf-8") as out_file:
        out_file.write("fid\tiid\tpop\tsupervised_label\toriginal_order\n")
        for row in sample_metadata_rows:
            out_file.write(
                f"{row['fid']}\t{row['iid']}\t{row['pop']}\t"
                f"{row['supervised_label']}\t{row['original_order']}\n"
            )

    # generate true ancestry from tree sequence using tspop
    log_msg("writing local and global ancestry tables")
    census_time = args.admixture_time + args.census_time_offset
    ancestry_table = get_local_ancestry_table(ts, census_time)
    sample_node_rows = []
    for ind_id in range(ts.num_individuals):
        nodes = sorted(int(node) for node in ts.individual(ind_id).nodes)
        for hap, node in enumerate(nodes, start=1):
            sample_node_rows.append((node, ind_id, hap))

    # for each population, generate pop specific ancestry
    for pop_idx, pop in enumerate(args.pops):
        start_node = 2 * pop_idx * args.sample_size
        end_node = 2 * ((pop_idx + 1) * args.sample_size) - 1
        pop_ancestry = ancestry_table[
            (ancestry_table["sample"] >= start_node)
            & (ancestry_table["sample"] <= end_node)
        ].copy()

        local_path = (
            f"{args.anc_dir}/{args.genetic_map}_{args.seed}_{pop}.tsv"
        )
        write_local_ancestry(
            pop_ancestry,
            sample_node_rows,
            local_path,
            args.chromosome,
            pop,
        )
        global_path = (
            f"{args.global_anc_dir}/{args.genetic_map}_{args.seed}_{pop}.tsv"
        )
        global_table = build_global_ancestry_table(
            pop_ancestry,
            sample_node_rows,
            pop,
        )
        write_global_ancestry(global_table, global_path)

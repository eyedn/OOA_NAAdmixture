###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           run_simulation.py
###############################################################################
# construct, simulate, and serialize OOA_NAAdmixture data products.


##### set up ##################################################################
import math
import pickle
import msprime
import stdpopsim
import tskit
from misc_utils import log_msg
from ancestry_utils.build_global_ancestry_table import (
    build_global_ancestry_table
)
from ancestry_utils.get_local_ancestry_table import (
    get_local_ancestry_table
)
from ancestry_utils.write_global_ancestry import (
    write_global_ancestry
)
from ancestry_utils.write_local_ancestry import (
    write_local_ancestry
)


##### internal functions ######################################################
'''
internal: parse a comma-delimited numeric CLI vector. Returns floating-point
values in the exact generation order supplied by the submitter.
'''
def _parse_float_list(raw_values):
    return [float(value) for value in raw_values.split(",") if value != ""]


'''
internal: generate metadata rows for supervised ADMIXTURE. Returns sample and
population IDs, supervised labels, and original sample order, using "-" for
ADX samples.
'''
def _build_metadata(pops, sample_size):
    rows = []
    original_order = 1
    for pop in pops:
        supervised_label = "-" if pop == "ADX" else pop
        for sample_idx in range(1, sample_size + 1):
            rows.append(
                {
                    "fid": "0",
                    "iid": f"{pop}_{sample_idx}",
                    "pop": pop,
                    "supervised_label": supervised_label,
                    "original_order": original_order,
                }
            )
            original_order += 1
    return rows


'''
internal: build and validate the OOA_NAAdmixture demography. Returns the
demography and the metadata required to reproduce its parameterization.
'''
def _build_demography(
    generation_time=25,
    mutation_rate=2.36e-8,
    t_af_years=148000,
    t_ooa_years=51000,
    t_eu0_years=23000,
    t_eg_years=5115,
    r_eu0=0.00307,
    r_eu=0.0195,
    r_af=0.0166,
    n_a=7310,
    n_af1=14474,
    n_b=1861,
    n_eu0=1032,
    m_af_b=15e-5,
    m_af_eu=2.5e-5,
    admixture_time=14.0,
    admix_generation_count=15,
    admix_mixing_generation_count=10,
    admix_ne_by_generation=None,
    admix_afr_props_by_generation=None,
    admix_eur_props_by_generation=None,
    admix_prioradmix_props_by_generation=None,
    admix_modern_growth_rate=0.023175,
    census_time_offset=1e-6,
):
    # define defaults for admixture process
    if admix_ne_by_generation is None:
        admix_ne_by_generation = [
            493.7874, 5755.8703, 15296.0328, 58666.1031,
            146967.1421, 260312.8139, 436921.0291, 858109.3961,
            1845355.3309, 3322796.2758, 3400701.3175, 3480432.8916,
            3562033.8225, 3645547.9384, 3731020.0950,
        ]
    if admix_afr_props_by_generation is None:
        admix_afr_props_by_generation = [
            0.850000, 0.904820, 0.791384, 0.786692, 0.719992,
            0.494960, 0.130456, 0.060000, 0.060000, 0.060000,
            0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
        ]
    if admix_eur_props_by_generation is None:
        admix_eur_props_by_generation = [
            0.150000, 0.080000, 0.080000, 0.080000, 0.080000,
            0.080000, 0.080000, 0.030000, 0.030000, 0.030000,
            0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
        ]
    if admix_prioradmix_props_by_generation is None:
        admix_prioradmix_props_by_generation = [
            0.000000, 0.015180, 0.128616, 0.133308, 0.200008,
            0.425040, 0.789544, 0.910000, 0.910000, 0.910000,
            1.000000, 1.000000, 1.000000, 1.000000, 1.000000,
        ]

    admix_vectors = [
        ("admix_ne_by_generation", admix_ne_by_generation),
        ("admix_afr_props_by_generation", admix_afr_props_by_generation),
        ("admix_eur_props_by_generation", admix_eur_props_by_generation),
        (
            "admix_prioradmix_props_by_generation",
            admix_prioradmix_props_by_generation,
        ),
    ]

    # validate admixture process inputs
    if admix_generation_count <= 0:
        raise ValueError("admix_generation_count must be positive")
    if admix_mixing_generation_count <= 0:
        raise ValueError("admix_mixing_generation_count must be positive")
    if admix_mixing_generation_count >= admix_generation_count:
        raise ValueError(
            "admix_mixing_generation_count must be less than "
            "admix_generation_count"
        )
    if admixture_time <= 0:
        raise ValueError("admixture_time must be positive")
    if census_time_offset <= 0:
        raise ValueError("census_time_offset must be positive")
    if not math.isclose(
        admixture_time,
        admix_generation_count - 1,
        rel_tol=0,
        abs_tol=1e-9,
    ):
        raise ValueError(
            "admixture_time must equal admix_generation_count - 1"
        )
    if admix_modern_growth_rate < 0:
        raise ValueError("admix_modern_growth_rate must be non-negative")
    for vector_name, vector_value in admix_vectors:
        if len(vector_value) != admix_generation_count:
            raise ValueError(
                f"{vector_name} must have length admix_generation_count"
            )
        if not all(math.isfinite(value) for value in vector_value):
            raise ValueError(f"{vector_name} must contain finite values")
    if not all(value > 0 for value in admix_ne_by_generation):
        raise ValueError("admix_ne_by_generation must contain positive values")
    for idx in range(admix_generation_count):
        props = [
            admix_afr_props_by_generation[idx],
            admix_eur_props_by_generation[idx],
            admix_prioradmix_props_by_generation[idx],
        ]
        if not all(value >= 0 for value in props):
            raise ValueError("admixture proportions must be non-negative")
        if not math.isclose(sum(props), 1.0, rel_tol=0, abs_tol=1e-6):
            raise ValueError(
                "admixture proportions must sum to 1 for every generation"
            )

    # times and population sizes from Tennessen 2012
    t_af = t_af_years / generation_time
    t_ooa = t_ooa_years / generation_time
    t_eu0 = t_eu0_years / generation_time
    t_eg = t_eg_years / generation_time
    n_eu1 = n_eu0 / math.exp(-r_eu0 * (t_eu0 - t_eg))
    n_eu = n_eu1 / math.exp(-r_eu * t_eg)
    n_af = n_af1 / math.exp(-r_af * t_eg)

    # initialize demography with Tennessen AFR/EUR and derived modern ADX.
    demography = msprime.Demography()
    demography.add_population(
        name="AFR",
        description="African population from the Tennessen 2012",
        initial_size=n_af,
        growth_rate=r_af,
        initially_active=True,
    )
    demography.add_population(
        name="EUR",
        description="European population from the Tennessen 2012",
        initial_size=n_eu,
        growth_rate=r_eu,
    )
    demography.add_population(
        name="ADX",
        description=(
            "Admixed population inspired from Hacker 2020 and Mooney 2023"
        ),
        initial_size=admix_ne_by_generation[-1],
        growth_rate=admix_modern_growth_rate,
    )

    # define intermediate ADX populations
    for gen_number in range(1, admix_mixing_generation_count):
        demography.add_population(
            name=f"ADX_G{gen_number}",
            description=(
                "One-generation Admixed population "
                f"for generation {gen_number}"
            ),
            initial_size=admix_ne_by_generation[gen_number - 1],
            growth_rate=0,
            initially_active=True,
        )

    # define the admixture event that derives modern ADX
    gen_number = admix_mixing_generation_count
    source_name = f"ADX_G{gen_number - 1}"
    demography.add_admixture(
        time=admix_generation_count - gen_number,
        derived="ADX",
        ancestral=["AFR", "EUR", source_name],
        proportions=[
            admix_afr_props_by_generation[gen_number - 1],
            admix_eur_props_by_generation[gen_number - 1],
            admix_prioradmix_props_by_generation[gen_number - 1],
        ],
    )

    # define intermediate admixture events for ADX_G* populations after ADX_G1.
    for gen_number in range(admix_mixing_generation_count - 1, 1, -1):
        raw_sources = ["AFR", "EUR", f"ADX_G{gen_number - 1}"]
        raw_proportions = [
            admix_afr_props_by_generation[gen_number - 1],
            admix_eur_props_by_generation[gen_number - 1],
            admix_prioradmix_props_by_generation[gen_number - 1],
        ]
        sources = [
            source for source, proportion in zip(raw_sources, raw_proportions)
            if proportion > 0
        ]
        proportions = [
            proportion for proportion in raw_proportions if proportion > 0
        ]
        demography.add_admixture(
            time=admix_generation_count - gen_number,
            derived=f"ADX_G{gen_number}",
            ancestral=sources,
            proportions=proportions,
        )

    # add a backward-time census after the ADX_G1 event for true local ancestry.
    demography.add_census(time=admixture_time + census_time_offset)
    demography.add_admixture(
        time=admixture_time,
        derived="ADX_G1",
        ancestral=["AFR", "EUR"],
        proportions=[
            admix_afr_props_by_generation[0],
            admix_eur_props_by_generation[0],
        ],
    )

    # define AFR/EUR OOA events from Tennessen 2012
    demography.set_symmetric_migration_rate(["AFR", "EUR"], m_af_eu)
    demography.add_symmetric_migration_rate_change(
        time=t_eg,
        populations=["AFR", "EUR"],
        rate=m_af_eu,
    )
    demography.add_population_parameters_change(
        time=t_eg,
        population="EUR",
        growth_rate=r_eu0,
        initial_size=n_eu1,
    )
    demography.add_population_parameters_change(
        time=t_eg,
        population="AFR",
        growth_rate=0,
        initial_size=n_af1,
    )
    demography.add_symmetric_migration_rate_change(
        time=t_eu0,
        populations=["AFR", "EUR"],
        rate=m_af_b,
    )
    demography.add_population_parameters_change(
        time=t_eu0,
        population="EUR",
        initial_size=n_b,
        growth_rate=0,
    )
    demography.add_population_split(
        time=t_ooa,
        derived=["EUR"],
        ancestral="AFR",
    )
    demography.add_migration_rate_change(time=t_ooa, rate=0)
    demography.add_population_parameters_change(
        time=t_af,
        population="AFR",
        initial_size=n_a,
    )

    # sort events and return the demography object and msprime metadata
    demography.sort_events()
    metadata = {
        "id": "OutOfAfrica_NorthAmericanAdmixture",
        "description": (
            "Tennessen 2012 with generation-indexed "
            "African American admixture"
        ),
        "generation_time": generation_time,
        "mutation_rate": mutation_rate,
        "admix_generation_count": admix_generation_count,
        "admix_mixing_generation_count": admix_mixing_generation_count,
        "admix_present_size": admix_ne_by_generation[-1],
        "admix_size_at_mixing_boundary": (
            admix_ne_by_generation[admix_mixing_generation_count - 1]
        ),
        "admix_modern_growth_rate": admix_modern_growth_rate,
    }
    return demography, metadata


##### main function ###########################################################
'''
simulate one OOA_NAAdmixture chromosome and write every downstream handoff:
tree sequence, demography metadata, VCF, sample metadata, and ancestry tables.
'''
def run_simulation(args):
    # define species and chromosome using stdpopsim
    species = stdpopsim.get_species("HomSap")
    contig = species.get_contig(
        f"chr{args.chromosome}", genetic_map=args.genetic_map
        )

    # define sample set for msprime simulator
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

    # build demography from above internal function
    demography, model_metadata = _build_demography(
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
        admix_ne_by_generation=_parse_float_list(
            args.admix_ne_by_generation
        ),
        admix_afr_props_by_generation=_parse_float_list(
            args.admix_afr_props_by_generation
        ),
        admix_eur_props_by_generation=_parse_float_list(
            args.admix_eur_props_by_generation
        ),
        admix_prioradmix_props_by_generation=_parse_float_list(
            args.admix_prioradmix_props_by_generation
        ),
        admix_modern_growth_rate=args.admix_modern_growth_rate,
        census_time_offset=args.census_time_offset,
    )

    # simulate the above demography with msprime and output the tree sequence
    log_msg(
        f"simulating tree sequence rep={args.seed} chr={args.chromosome}"
    )

    # generate tree
    ts = msprime.sim_ancestry(
        samples=sample_sets,
        demography=demography,
        recombination_rate=contig.recombination_map,
        sequence_length=contig.length,
        ploidy=2,
        random_seed=args.seed,
        model=args.msprime_model,
    )

    # generate mutations
    ts = msprime.sim_mutations(
        ts,
        rate=args.mutation_rate,
        random_seed=args.seed + 1,
    )

    # verify the result is a valid tree sequence
    if not isinstance(ts, tskit.TreeSequence):
        raise TypeError("Expected a tskit.TreeSequence result")
    ts_path = f"{args.tree_prefix}.ts"
    ts.dump(ts_path)

    # pickle demography and demography metadata objects
    log_msg(
        f"writing demography and demography metadata pickles rep={args.seed} "
        f"chr={args.chromosome}"
    )
    with open(f"{args.pickle_prefix}.demography.pkl", "wb") as out_file:
        pickle.dump(demography, out_file)
    with open(f"{args.pickle_prefix}.model_metadata.pkl", "wb") as out_file:
        pickle.dump(model_metadata, out_file)

    # output vcf resulting from tree sequence
    log_msg(
        f"writing all-sample VCF rep={args.seed} chr={args.chromosome} "
        f"path={args.vcf_path}"
    )
    with open(args.vcf_path, "w", encoding="utf-8") as out_file:
        ts.write_vcf(
            out_file,
            contig_id=f"chr{args.chromosome}",
            individual_names=sample_names,
            position_transform=lambda positions: [pos + 1 for pos in positions],
        )

    # output sample metadata for downstream "admixture --supervised"
    log_msg(
        f"writing sample metadata rep={args.seed} chr={args.chromosome} "
        f"path={args.sample_metadata_path}"
    )
    sample_metadata_rows = _build_metadata(
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
    log_msg(
        f"writing local and global ancestry tables rep={args.seed} "
        f"chr={args.chromosome}"
    )
    census_time = args.admixture_time + args.census_time_offset
    ancestry_table = get_local_ancestry_table(ts, census_time)
    sample_node_rows = []
    for ind_id in range(ts.num_individuals):
        nodes = sorted(int(node) for node in ts.individual(ind_id).nodes)
        for hap, node in enumerate(nodes, start=1):
            sample_node_rows.append((node, ind_id, hap))

    # write local and global tspop ancestry outputs for each population.
    for pop_idx, pop in enumerate(args.pops):
        start_node = 2 * pop_idx * args.sample_size
        end_node = 2 * ((pop_idx + 1) * args.sample_size) - 1
        pop_ancestry = ancestry_table[
            (ancestry_table["sample"] >= start_node)
            & (ancestry_table["sample"] <= end_node)
        ].copy()

        local_path = (
            f"{args.anc_dir}/{args.genetic_map}_{args.seed}_"
            f"chr{args.chromosome}_{pop}.tsv"
        )
        write_local_ancestry(
            pop_ancestry,
            sample_node_rows,
            local_path,
            args.chromosome,
            pop,
        )
        global_path = (
            f"{args.global_anc_dir}/{args.genetic_map}_{args.seed}_"
            f"chr{args.chromosome}_{pop}.tsv"
        )
        global_table = build_global_ancestry_table(
            pop_ancestry,
            sample_node_rows,
            pop,
        )
        write_global_ancestry(global_table, global_path)

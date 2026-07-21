###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           calc_stats.py
###############################################################################

# overview: calculate simulation stats. from tree-sequence and PLINK inputs.


##### set up ##################################################################
from itertools import combinations
from pathlib import Path
import csv
import allel
import numpy as np
import pandas as pd
import tszip
from misc_utils import log_msg
from .parse_king_file import parse_king_file
from .read_fam_order import read_q_rows


##### internal functions ######################################################
ANCESTRY_COLUMNS = [
    "rep",
    "chrom",
    "pop",
    "sample_id",
    "vcf_sample_id",
    "afr_tspop",
    "eur_tspop",
    "afr_q",
    "eur_q",
    "span",
]

'''
internal: return sample nodes for one population in tree-sequence sample order.
'''
def _sample_nodes_for_pop(ts, pops, sample_size, pop):
    pop_idx = pops.index(pop)
    start = pop_idx * sample_size
    end = (pop_idx + 1) * sample_size
    nodes = []
    for ind_id in range(start, end):
        nodes.extend(int(node) for node in ts.individual(ind_id).nodes)
    return nodes


'''
internal: build one-dimensional SFS rows for every population. Returns derived
allele-count bins for one chromosome and replicate.
'''
def _build_1d_sfs_rows(ts, rep, pops, sample_size, chrom=None):
    rows = []
    for pop in pops:
        sample_nodes = _sample_nodes_for_pop(ts, pops, sample_size, pop)
        spectrum = ts.allele_frequency_spectrum(
            [sample_nodes],
            polarised=True,
            span_normalise=False,
            mode="site",
        )
        for derived_count, count in enumerate(spectrum):
            rows.append(
                {
                    "rep": rep,
                    "chrom": chrom,
                    "pop": pop,
                    "derived_allele_count": derived_count,
                    "count": count,
                }
            )
    return rows


'''
internal: build pairwise two-dimensional SFS rows. Returns joint derived-count
bins for every population pair.
'''
def _build_2d_sfs_rows(ts, rep, pops, sample_size, chrom=None):
    rows = []
    for pop1, pop2 in combinations(pops, 2):
        sample_sets = [
            _sample_nodes_for_pop(ts, pops, sample_size, pop1),
            _sample_nodes_for_pop(ts, pops, sample_size, pop2),
        ]
        spectrum = ts.allele_frequency_spectrum(
            sample_sets,
            polarised=True,
            span_normalise=False,
            mode="site",
        )
        for pop1_count in range(spectrum.shape[0]):
            for pop2_count in range(spectrum.shape[1]):
                rows.append(
                    {
                        "rep": rep,
                        "chrom": chrom,
                        "pop1": pop1,
                        "pop2": pop2,
                        "pop1_count": pop1_count,
                        "pop2_count": pop2_count,
                        "count": spectrum[pop1_count, pop2_count],
                    }
                )
    return rows


'''
internal: build separate pi and Watterson theta rows for each population using
the statistic-specific estimators and one explicit sequence span.
'''
def _build_pi_theta_rows(ts, rep, pops, sample_size, mutation_rate, chrom=None):
    rows = []
    for pop in pops:
        sample_nodes = _sample_nodes_for_pop(ts, pops, sample_size, pop)

        # pi calculation
        pi_value = ts.diversity([sample_nodes], mode="site")[0]

        # theta calculation
        wattersons_const = sum(
            1.0 / value for value in range(1, len(sample_nodes))
        )
        segregating_sites = ts.segregating_sites(
            [sample_nodes],
            mode="site",
        )[0]
        theta_value = segregating_sites / wattersons_const

        # append statistics to returned data rows
        rows.append(
            {
                "rep": rep,
                "chrom": chrom,
                "pop": pop,
                "stat": "pi",
                "value": pi_value,
                "ne_value": pi_value / (4 * mutation_rate),
                "mutation_rate": mutation_rate,
                "span": float(ts.sequence_length),
                "segregating_sites": None,
                "wattersons_const": wattersons_const,
            }
        )
        rows.append(
            {
                "rep": rep,
                "chrom": chrom,
                "pop": pop,
                "stat": "theta",
                "value": theta_value,
                "ne_value": theta_value / (4 * mutation_rate),
                "mutation_rate": mutation_rate,
                "span": float(ts.sequence_length),
                "segregating_sites": segregating_sites,
                "wattersons_const": wattersons_const,
            }
        )
    return rows


'''
internal: build LD-decay summary rows within non-overlapping genome windows:
    - convert tree-sequence genotypes to diploid genotype arrays
    - filter variants by minor allele frequency
    - divide the chromosome into non-overlapping windows
    - calculate pairwise Rogers–Huff correlations with scikit-allel
    - square them to obtain r2
    - bin each r2 values by physical distance of the respective loci pair
    - record mean and summed r2 plus pair counts
'''
def _build_ld_decay_rows(
    ts,
    rep,
    pops,
    sample_size,
    chrom=None,
    window_size_bp=2_000_000,
    distance_bin_bp=5_000,
    maf_threshold=0.10,
):
    rows = []
    positions = np.array([site.position for site in ts.sites()])

    # generate the haploid genotype matrix from the tree sequence.
    haploid_genotypes = ts.genotype_matrix().astype(int)
    if haploid_genotypes.shape[1] % 2 != 0:
        raise ValueError("Expected diploid genotypes with an even node count")

    # convert haplotype genotype matrix to diploid genotype matrix
    num_individuals = haploid_genotypes.shape[1] // 2
    genotypes = haploid_genotypes.reshape(
        haploid_genotypes.shape[0],
        num_individuals,
        2,
    )
    requested_individuals = len(pops) * sample_size
    if requested_individuals > num_individuals:
        raise ValueError(
            "Requested populations require more individuals than available"
        )

    # convert matrix to allel.GenotypeArray object
    genotype_array = allel.GenotypeArray(genotypes[:, :requested_individuals])
    allele_counts = genotype_array.count_alleles()
    alt_freq = allele_counts[:, 1] / allele_counts.sum(axis=1)
    maf = np.minimum(alt_freq, 1 - alt_freq)
    keep = maf >= maf_threshold
    positions = positions[keep]

    # generate alternative allele counts of loci that passed the MAF threshold
    alt_counts = genotype_array[keep].to_n_alt()

    # calculate LD decay separately for each population's sample columns.
    for pop in pops:
        # define individual indexes for the target population.
        pop_idx = pops.index(pop)
        start = pop_idx * sample_size
        end = start + sample_size
        pop_alt_counts = alt_counts[:, start:end]

        # calculate all pairwise r2 within each non-overlapping window
        for window_start in range(0, int(ts.sequence_length), window_size_bp):
            window_end = min(window_start + window_size_bp, ts.sequence_length)
            in_window = np.where(
                (positions >= window_start) & (positions < window_end)
            )[0]
            if len(in_window) < 2:
                continue

            # define positions in the window
            window_positions = positions[in_window]
            block = pop_alt_counts[in_window, :]

            # calculate r2 and retain pairs that fall within the window.
            r_values = allel.rogers_huff_r(block)
            row_idx, col_idx = np.triu_indices(len(window_positions), k=1)
            distances = window_positions[col_idx] - window_positions[row_idx]
            r2_values = np.asarray(r_values) ** 2
            valid = (distances > 0) & (distances <= window_size_bp)
            distances = distances[valid]
            r2_values = r2_values[valid]

            # bin r2 values by the distance of the respective loci pair
            bin_values = {}
            distance_bins = (
                ((distances.astype(int) - 1) // distance_bin_bp) + 1
            ) * distance_bin_bp
            for distance_bin, r2_value in zip(distance_bins, r2_values):
                bin_values.setdefault(int(distance_bin), []).append(
                    float(r2_value)
                )

            # summarize r2 values by window and bin
            for distance_bin, r2_values in sorted(bin_values.items()):
                rows.append(
                    {
                        "rep": rep,
                        "chrom": chrom,
                        "pop": pop,
                        "window_start": int(window_start),
                        "window_end": int(window_end),
                        "distance_bin_bp": int(distance_bin),
                        "mean_r2": float(np.nanmean(r2_values)),
                        "sum_r2": float(np.nansum(r2_values)),
                        "n_pairs": len(r2_values),
                    }
                )
    return rows


'''
internal: join tspop ancestry and supervised ADMIXTURE estimates by sample.
Returns rows with both tree-sequence and VCF identifiers.
'''
def _build_ancestry_rows(
    rep,
    pops,
    genetic_map,
    global_anc_dir,
    q_path,
    fam_path,
    chrom=None,
):
    global_rows = []
    for pop in pops:
        global_path = Path(global_anc_dir) / (
            f"{genetic_map}_{rep}_chr{chrom}_{pop}.tsv"
        )
        with open(global_path, "r", encoding="utf-8") as in_file:
            reader = csv.DictReader(in_file, delimiter="\t")
            for row in reader:
                global_rows.append(
                    {
                        "rep": rep,
                        "chrom": chrom,
                        "pop": pop,
                        "sample_id": row["sample_id"],
                        "vcf_sample_id": row["vcf_sample_id"],
                        "afr_tspop": float(row["AFR_prop"]),
                        "eur_tspop": float(row["EUR_prop"]),
                        "span": float(row.get("span", 1.0)),
                    }
                )

    q_by_sample = {
        (row["pop"], row["vcf_sample_id"]): row
        for row in read_q_rows(rep, q_path, fam_path)
    }
    ancestry_rows = []
    for row in global_rows:
        q_row = q_by_sample.get((row["pop"], row["vcf_sample_id"]), {})
        ancestry_rows.append(
            {
                **row,
                "afr_q": q_row.get("afr_q"),
                "eur_q": q_row.get("eur_q"),
            }
        )

    return ancestry_rows


'''
internal: build a pandas table with the canonical joined-ancestry columns.
'''
def _build_ancestry_table(*args, **kwargs):
    rows = _build_ancestry_rows(*args, **kwargs)
    return pd.DataFrame(rows, columns=ANCESTRY_COLUMNS)


##### main function ###########################################################
'''
calculate all simulation statistics for one chromosome and replicate. Writes
joined ancestry, population and combined KING tables, pi, theta, one- and
two-dimensional SFS tables, and LD-decay summaries.
'''
def calc_stats(args):
    stats_dir = Path(args.stats_dir)
    stats_dir.mkdir(parents=True, exist_ok=True)

    # load compressed tree sequence
    tree_tsz_path = Path(args.tree_tsz_path)
    if not str(tree_tsz_path).endswith(".ts.tsz"):
        raise ValueError(
            f"Tree sequence input must end with .ts.tsz: {tree_tsz_path}"
        )
    if not tree_tsz_path.exists():
        raise FileNotFoundError(
            f"Missing compressed tree sequence {tree_tsz_path}"
        )
    log_msg(
        f"loading tree sequence rep={args.rep} chr={args.chr} "
        f"path={tree_tsz_path}"
    )
    ts = tszip.decompress(tree_tsz_path)
    log_msg(f"loaded tree sequence rep={args.rep} chr={args.chr}")

    # join tspop and supervised ADMIXTURE ancestry results.
    q_path = Path(args.admixture_dir) / (
        f"{args.genetic_map}_{args.rep}_chr{args.chr}_all.2.Q"
    )
    log_msg(f"parsing ADMIXTURE Q rep={args.rep} chr={args.chr} path={q_path}")
    ancestry_table = _build_ancestry_table(
        rep=args.rep,
        pops=args.pops,
        genetic_map=args.genetic_map,
        global_anc_dir=args.global_anc_dir,
        q_path=q_path,
        fam_path=args.admixture_fam_path,
        chrom=args.chr,
    )
    ancestry_table.to_csv(
        stats_dir / f"ancestry.rep_{args.rep}.chr{args.chr}.tsv",
        sep="\t",
        index=False
    )
    ancestry_table.to_parquet(
        stats_dir / f"ancestry.rep_{args.rep}.chr{args.chr}.parquet",
        index=False
    )

    # standardize population KING outputs and write one combined kinship table.
    king_tables = []
    for pop in args.pops:
        king_path = Path(args.king_dir) / (
            f"{args.genetic_map}_{args.rep}_chr{args.chr}_{pop}.kin0"
        )
        king_table = parse_king_file(king_path, args.rep, pop)
        king_table.insert(1, "chrom", args.chr)
        king_table.to_csv(
            Path(args.king_dir) / (
                f"{args.genetic_map}_{args.rep}_chr{args.chr}_{pop}.tsv"
            ),
            sep="\t",
            index=False
        )
        king_tables.append(king_table)
    king_combined = pd.concat(king_tables, ignore_index=True)
    king_combined.to_csv(
        stats_dir / f"kinship.rep_{args.rep}.chr{args.chr}.tsv",
        sep="\t",
        index=False
    )
    king_combined.to_parquet(
        stats_dir / f"kinship.rep_{args.rep}.chr{args.chr}.parquet",
        index=False
    )

    # derive and process pi and theta from the tree sequence
    log_msg(f"calculating pi/theta rep={args.rep} chr={args.chr}")
    pi_theta = pd.DataFrame(
        _build_pi_theta_rows(
            ts,
            args.rep,
            args.pops,
            args.sample_size,
            args.mutation_rate,
            args.chr,
        )
    )

    # derive and process 1D and 2D sfs from the tree sequence
    log_msg(f"calculating SFS rep={args.rep} chr={args.chr}")
    sfs = pd.DataFrame(
        _build_1d_sfs_rows(
            ts,
            args.rep,
            args.pops,
            args.sample_size,
            args.chr,
        )
    )
    sfs_2d = pd.DataFrame(
        _build_2d_sfs_rows(
            ts,
            args.rep,
            args.pops,
            args.sample_size,
            args.chr,
        )
    )

    # derive and process ld decay from the tree sequence
    log_msg(f"calculating LD decay rep={args.rep} chr={args.chr}")
    ld_decay = pd.DataFrame(
        _build_ld_decay_rows(
            ts,
            args.rep,
            args.pops,
            args.sample_size,
            args.chr,
            args.ld_decay_window_size_bp,
            args.ld_decay_distance_bin_bp,
            args.ld_decay_maf_threshold,
        )
    )

    # combine summary tables into one chromosome statistics table.
    tables = {
        "pi_theta_stats": pi_theta,
        "sfs": sfs,
        "sfs_2d": sfs_2d,
        "ld_decay": ld_decay,
    }
    for table_name, table in tables.items():
        table.to_csv(
            stats_dir / f"{table_name}.rep_{args.rep}.chr{args.chr}.tsv",
            sep="\t",
            index=False,
        )
        table.to_parquet(
            stats_dir / f"{table_name}.rep_{args.rep}.chr{args.chr}.parquet",
            index=False,
        )
    log_msg(f"wrote chromosome stats rep={args.rep} chr={args.chr}")

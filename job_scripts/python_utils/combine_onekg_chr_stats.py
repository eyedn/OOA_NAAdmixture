###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           combine_onekg_chr_stats.py
###############################################################################

# overview: combine 1000 Genomes chromosome stats. across population outputs.


##### set up ##################################################################
from itertools import combinations
from pathlib import Path
import argparse
from onekg_utils.build_folded_2d_sfs_rows import build_folded_2d_sfs_rows
from onekg_utils.read_tsv_rows import read_tsv_rows
from onekg_utils.write_stats_table import write_stats_table


##### arguments ###############################################################
'''
define the statistics directory, target chromosome, and ordered population
labels used to combine per-population empirical outputs.
'''
parser = argparse.ArgumentParser()
parser.add_argument("--stats-dir", required=True)
parser.add_argument("--chrom", required=True)
parser.add_argument("--pops", nargs="+", required=True)


##### main ####################################################################
if __name__ == "__main__":
    args = parser.parse_args()
    stats_dir = Path(args.stats_dir)
    # combine population-specific tables that already share one schema.
    partial_tables = [
        "pi_theta_stats_intergenic",
        "pi_theta_stats_full_callable_chrom",
        "sfs",
        "ld_decay",
        "kinship_unrelated",
    ]
    combined = {}
    for table_name in partial_tables:
        rows = []
        for pop in args.pops:
            path = stats_dir / (
                f"{table_name}.rep_0.chr{args.chrom}.{pop}.tsv"
            )
            rows.extend(read_tsv_rows(path))
        combined[table_name] = rows
        write_stats_table(
            stats_dir / f"{table_name}.rep_0.chr{args.chrom}",
            rows,
        )
        write_stats_table(
            stats_dir / f"{table_name}.chr{args.chrom}",
            rows,
        )
    write_stats_table(
        stats_dir / f"kinship.rep_0.chr{args.chrom}",
        combined["kinship_unrelated"],
    )
    write_stats_table(
        stats_dir / f"kinship.chr{args.chrom}",
        combined["kinship_unrelated"],
    )

    # build the folded pairwise SFS from shared per-population allele counts.
    counts_by_pop = {}
    for pop in args.pops:
        path = stats_dir / (
            f"allele_counts.rep_0.chr{args.chrom}.{pop}.tsv"
        )
        counts_by_pop[pop] = {
            int(row["position"]): (
                int(row["ref_count"]),
                int(row["alt_count"]),
            )
            for row in read_tsv_rows(path)
        }
    sfs_2d = []
    for pop1, pop2 in combinations(args.pops, 2):
        sfs_2d.extend(
            build_folded_2d_sfs_rows(
                0,
                args.chrom,
                pop1,
                pop2,
                counts_by_pop[pop1],
                counts_by_pop[pop2],
            )
        )
    write_stats_table(
        stats_dir / f"sfs_2d.rep_0.chr{args.chrom}",
        sfs_2d,
    )
    write_stats_table(
        stats_dir / f"sfs_2d.chr{args.chrom}",
        sfs_2d,
    )

    # write chromosome-level ancestry and variant-QC handoffs.
    ancestry_tables = [
        "ancestry_ADMIXTURE_super",
        "ancestry_ADMIXTURE_multik",
        "ancestry_fastStructure_multik",
        "fastStructure_chooseK",
    ]
    for table_name in ancestry_tables:
        rows = read_tsv_rows(
            stats_dir / f"{table_name}.rep_0.chr{args.chrom}.tsv"
        )
        write_stats_table(
            stats_dir / f"{table_name}.chr{args.chrom}",
            rows,
        )
        write_stats_table(
            stats_dir / f"{table_name}.rep_0.chr{args.chrom}",
            rows,
        )
    for old_family in ("ancestry", "ancestry_multik"):
        old_prefix = stats_dir / f"{old_family}.chr{args.chrom}"
        for extension in ("tsv", "parquet"):
            Path(f"{old_prefix}.{extension}").unlink(missing_ok=True)

    qc_rows = read_tsv_rows(
        stats_dir
        / f"variant_qc.rep_0.chr{args.chrom}.{args.pops[0]}.tsv"
    )
    qc_row = {**qc_rows[0], "pop": "ALL"}
    write_stats_table(
        stats_dir / f"variant_qc.chr{args.chrom}",
        [qc_row],
    )
    for key, value in qc_row.items():
        if key not in {"rep", "chrom", "pop"}:
            print(f"variant_qc chr={args.chrom} {key}={value}")

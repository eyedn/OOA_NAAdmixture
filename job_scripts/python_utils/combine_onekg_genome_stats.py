###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           combine_onekg_genome_stats.py
###############################################################################

# overview: combine 1000 Genomes genome stats. across population outputs.


##### set up ##################################################################
from pathlib import Path
import argparse
from onekg_utils.aggregate_folded_2d_sfs_rows import (
    aggregate_folded_2d_sfs_rows,
)
from onekg_utils.aggregate_variant_qc_rows import aggregate_variant_qc_rows
from onekg_utils.read_tsv_rows import read_tsv_rows
from onekg_utils.write_stats_table import write_stats_table


##### arguments ###############################################################
'''
define the statistics directory and ordered chromosome and population lists
used to combine empirical genome outputs.
'''
parser = argparse.ArgumentParser()
parser.add_argument("--stats-dir", required=True)
parser.add_argument("--chroms", nargs="+", required=True)
parser.add_argument("--pops", nargs="+", required=True)


##### main ####################################################################
if __name__ == "__main__":
    args = parser.parse_args()
    stats_dir = Path(args.stats_dir)
    # combine population-specific genome tables that share one schema.
    pop_tables = [
        "pi_theta_stats_intergenic",
        "pi_theta_stats_full_callable_chrom",
        "sfs",
        "ld_decay",
        "kinship_unrelated",
    ]
    combined = {}
    for table_name in pop_tables:
        rows = []
        for pop in args.pops:
            rows.extend(
                read_tsv_rows(stats_dir / f"{table_name}.rep_0.{pop}.tsv")
            )
        combined[table_name] = rows
        write_stats_table(stats_dir / f"{table_name}.rep_0", rows)
        write_stats_table(stats_dir / table_name, rows)
    write_stats_table(
        stats_dir / "kinship.rep_0",
        combined["kinship_unrelated"],
    )
    write_stats_table(
        stats_dir / "kinship",
        combined["kinship_unrelated"],
    )

    # aggregate pairwise SFS and variant-QC rows across chromosomes.
    sfs_2d_rows = []
    qc_rows = []
    for chrom in args.chroms:
        sfs_2d_rows.extend(
            read_tsv_rows(stats_dir / f"sfs_2d.rep_0.chr{chrom}.tsv")
        )
        qc_rows.extend(
            read_tsv_rows(stats_dir / f"variant_qc.chr{chrom}.tsv")
        )
    sfs_2d = aggregate_folded_2d_sfs_rows(sfs_2d_rows)
    variant_qc = aggregate_variant_qc_rows(qc_rows)
    for output_name, rows in (
        ("sfs_2d.rep_0", sfs_2d),
        ("sfs_2d", sfs_2d),
        ("variant_qc.rep_0", variant_qc),
        ("variant_qc", variant_qc),
    ):
        write_stats_table(stats_dir / output_name, rows)

    ancestry_tables = [
        "ancestry_ADMIXTURE_super",
        "ancestry_ADMIXTURE_multik",
        "ancestry_fastStructure_multik",
        "fastStructure_chooseK",
    ]
    for table_name in ancestry_tables:
        rows = read_tsv_rows(stats_dir / f"{table_name}.rep_0.tsv")
        write_stats_table(stats_dir / table_name, rows)
    for old_family in ("ancestry", "ancestry_multik"):
        for extension in ("tsv", "parquet"):
            (stats_dir / f"{old_family}.{extension}").unlink(
                missing_ok=True
            )

###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           write_stats_table.py
###############################################################################

# overview: write empirical statistics tables in canonical TSV and Parquet
# formats.

# pattern: Imperative Shell


##### set up ##################################################################
from pathlib import Path
import csv
import subprocess
import sys


EMPTY_TABLE_COLUMNS = {
    "allele_counts": [
        "rep", "chrom", "position", "pop", "ref_count", "alt_count",
    ],
    "ancestry": [
        "rep", "chrom", "pop", "sample_id", "vcf_sample_id",
        "afr_tspop", "eur_tspop", "afr_q", "eur_q", "span",
    ],
    "ancestry_multik": [
        "rep", "chrom", "pop", "sample_id", "vcf_sample_id", "k",
        "component_1_q", "component_2_q", "component_3_q",
        "component_4_q", "component_5_q", "span",
    ],
    "kinship": ["rep", "chrom", "pop", "id1", "id2", "kinship"],
    "kinship_unrelated": [
        "rep", "chrom", "pop", "id1", "id2", "kinship",
    ],
    "ld_decay": [
        "rep", "chrom", "pop", "window_start", "window_end",
        "distance_bin_bp", "mean_r2", "sum_r2", "n_pairs",
    ],
    "pi_theta_stats_intergenic": [
        "rep", "chrom", "pop", "stat", "value", "ne_value",
        "mutation_rate", "span", "segregating_sites",
        "wattersons_const",
    ],
    "pi_theta_stats_full_callable_chrom": [
        "rep", "chrom", "pop", "stat", "value", "ne_value",
        "mutation_rate", "span", "segregating_sites",
        "wattersons_const",
    ],
    "sfs": ["rep", "chrom", "pop", "minor_allele_count", "count"],
    "sfs_2d": [
        "rep", "chrom", "pop1", "pop2", "pop1_minor_allele_count",
        "pop2_minor_allele_count", "count",
    ],
}


##### main function ###########################################################
'''
write one statistics table as canonical TSV and Parquet files. Empty tables use
the registered schema so downstream combination retains stable columns.
'''
def write_stats_table(output_prefix, rows):
    prefix = Path(output_prefix)
    tsv_path = Path(f"{prefix}.tsv")
    parquet_path = Path(f"{prefix}.parquet")
    table_name = prefix.name.split(".", 1)[0]
    fieldnames = (
        list(rows[0])
        if rows
        else EMPTY_TABLE_COLUMNS.get(table_name, [])
    )
    if not fieldnames:
        raise ValueError(f"Cannot infer empty table schema for {prefix}")
    with open(tsv_path, "w", encoding="utf-8", newline="") as out_file:
        writer = csv.DictWriter(
            out_file,
            fieldnames=fieldnames,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)
    converter = Path(__file__).parents[1] / "write_parquet.py"
    subprocess.run(
        [sys.executable, str(converter), str(tsv_path), str(parquet_path)],
        check=True,
    )

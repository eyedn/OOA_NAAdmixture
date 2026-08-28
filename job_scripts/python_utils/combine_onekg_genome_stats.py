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
from onekg_utils.read_tsv_rows import read_tsv_rows
from shared_utils.write_stats_table import write_stats_table


##### internal functions ######################################################
'''
internal: sum chromosome variant-QC counts while retaining sample counts.
'''
def _aggregate_variant_qc_rows(rows):
    output = {"rep": 0, "chrom": "genome", "pop": "ALL"}
    for row in rows:
        for key, value in row.items():
            if key in {"rep", "chrom", "pop"} or value in {None, ""}:
                continue
            number = int(float(value))
            if key.startswith("retained_") and key.endswith("_samples"):
                output[key] = max(output.get(key, 0), number)
            else:
                output[key] = output.get(key, 0) + number
    return [output]


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
        "kinship_unrelated"
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
        combined["kinship_unrelated"]
    )
    write_stats_table(
        stats_dir / "kinship",
        combined["kinship_unrelated"]
    )

    # aggregate variant-QC rows across chromosomes.
    qc_rows = []
    for chrom in args.chroms:
        qc_rows.extend(
            read_tsv_rows(stats_dir / f"variant_qc.chr{chrom}.tsv")
        )
    variant_qc = _aggregate_variant_qc_rows(qc_rows)
    for output_name, rows in (
        ("variant_qc.rep_0", variant_qc),
        ("variant_qc", variant_qc)
    ):
        write_stats_table(stats_dir / output_name, rows)

    ancestry_tables = [
        "ancestry_ADMIXTURE_super",
        "ancestry_ADMIXTURE_multik",
        "ancestry_fastStructure_multik",
        "fastStructure_chooseK"
    ]
    for table_name in ancestry_tables:
        rows = read_tsv_rows(stats_dir / f"{table_name}.rep_0.tsv")
        write_stats_table(stats_dir / table_name, rows)
    for old_family in ("ancestry", "ancestry_multik"):
        for extension in ("tsv", "parquet"):
            (stats_dir / f"{old_family}.{extension}").unlink(
                missing_ok=True
            )

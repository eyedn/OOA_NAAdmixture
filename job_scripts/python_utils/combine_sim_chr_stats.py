###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           combine_sim_chr_stats.py
###############################################################################

# overview: combine simulated chromosome statistics across replicate outputs.


##### set up ##################################################################
from pathlib import Path
import argparse
import csv
import subprocess
import sys

'''
define supported chrom-level table types:
    - ancestry
    - kinship
    - unrelated kinship
    - pi/theta
    - one-dimensional SFS
    - LD decay
'''
TABLE_NAMES = [
    "ancestry",
    "ancestry_ADMIXTURE_multik",
    "ancestry_fastStructure_multik",
    "fastStructure_chooseK",
    "kinship",
    "kinship_unrelated",
    "pi_theta_stats",
    "sfs",
    "ld_decay",
]


##### internal functions #####################################################
''' internal: read and concatenate required TSV paths without skipping
input. '''
def _combine_table_paths(paths):
    missing_paths = [Path(path) for path in paths if not Path(path).exists()]
    if missing_paths:
        raise FileNotFoundError(", ".join(str(path) for path in missing_paths))
    rows = []
    for path in paths:
        with open(path, "r", encoding="utf-8", newline="") as in_file:
            rows.extend(csv.DictReader(in_file, delimiter="\t"))
    return rows


''' internal: write a combined TSV and invoke the canonical Parquet writer. '''
def _write_combined_stats_table(stats_dir, output_name, rows):
    stats_path = Path(stats_dir)
    tsv_path = stats_path / f"{output_name}.tsv"
    parquet_path = stats_path / f"{output_name}.parquet"
    fieldnames = list(rows[0]) if rows else []
    with open(tsv_path, "w", encoding="utf-8", newline="") as out_file:
        writer = csv.DictWriter(
            out_file,
            fieldnames=fieldnames,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)
    converter = Path(__file__).parent / "write_parquet.py"
    subprocess.run(
        [sys.executable, str(converter), str(tsv_path), str(parquet_path)],
        check=True,
    )


##### arguments ###############################################################
'''
define command-line arguments for
    - statistics output directory
    - number of replicates
    - list of chromosomes
    - chromosome index
    - list of tables to include / optional table selection
'''
parser = argparse.ArgumentParser()
parser.add_argument("--stats-dir", required=True)
parser.add_argument("--num-reps", type=int, required=True)
parser.add_argument("--chroms", nargs="+", required=True)
parser.add_argument("--chrom-index", type=int, required=True)
parser.add_argument(
    "--tables",
    nargs="+",
    choices=TABLE_NAMES,
    default=TABLE_NAMES
)


##### main ###################################################################
if __name__ == "__main__":
    args = parser.parse_args()
    if args.chrom_index < 1 or args.chrom_index > len(args.chroms):
        raise ValueError(
            f"--chrom-index must be between 1 and {len(args.chroms)}"
        )

    stats_path = Path(args.stats_dir)
    chrom = args.chroms[args.chrom_index - 1]
    for table_name in args.tables:
        paths = [
            stats_path / f"{table_name}.rep_{rep}.chr{chrom}.tsv"
            for rep in range(1, args.num_reps + 1)
        ]
        try:
            rows = _combine_table_paths(paths)
        except FileNotFoundError as exc:
            raise FileNotFoundError(
                f"Missing chromosome files for {table_name}: {exc}"
            ) from exc
        _write_combined_stats_table(
            stats_path,
            f"{table_name}.chr{chrom}",
            rows
        )

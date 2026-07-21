###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           combine_sim_chr_stats.py
###############################################################################
# combine one simulated chromosome statistic across replicate outputs.


##### set up ##################################################################
from pathlib import Path
import argparse
from stats_utils.combine_table_paths import combine_table_paths
from stats_utils.write_combined_stats_table import write_combined_stats_table

'''
define supported chrom-level table types:
    - ancestry
    - kinship
    - unrelated kinship
    - pi/theta
    - one-dimensional SFS
    - two-dimensional SFS
    - LD decay
'''
TABLE_NAMES = [
    "ancestry",
    "kinship",
    "kinship_unrelated",
    "pi_theta_stats",
    "sfs",
    "sfs_2d",
    "ld_decay",
]


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
    default=TABLE_NAMES,
)


##### main ####################################################################
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
            rows = combine_table_paths(paths)
        except FileNotFoundError as exc:
            raise FileNotFoundError(
                f"Missing chromosome files for {table_name}: {exc}"
            ) from exc
        write_combined_stats_table(
            stats_path,
            f"{table_name}.chr{chrom}",
            rows,
        )

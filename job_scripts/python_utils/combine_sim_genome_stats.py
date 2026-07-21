###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           combine_sim_genome_stats.py
###############################################################################

# overview: combine simulated genome statistics across replicate outputs.


##### set up ##################################################################
from pathlib import Path
import argparse
from sim_utils.combine_table_paths import combine_table_paths
from sim_utils.write_combined_stats_table import write_combined_stats_table


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
define the statistics directory, replicate count, and optional subset of
genome-level tables to combine.
'''
parser = argparse.ArgumentParser()
parser.add_argument("--stats-dir", required=True)
parser.add_argument("--num-reps", type=int, required=True)
parser.add_argument(
    "--tables",
    nargs="+",
    choices=TABLE_NAMES,
    default=TABLE_NAMES,
)


##### main ####################################################################
if __name__ == "__main__":
    args = parser.parse_args()
    stats_path = Path(args.stats_dir)
    for table_name in args.tables:
        paths = [
            stats_path / f"{table_name}.rep_{rep}.tsv"
            for rep in range(1, args.num_reps + 1)
        ]
        try:
            rows = combine_table_paths(paths)
        except FileNotFoundError as exc:
            raise FileNotFoundError(
                f"Missing per-replicate files for {table_name}: {exc}"
            ) from exc
        write_combined_stats_table(stats_path, table_name, rows)

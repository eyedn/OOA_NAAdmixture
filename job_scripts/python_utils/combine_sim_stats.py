###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           combine_sim_stats.py
###############################################################################


import argparse
from stats_utils.combine_sim_stats_tables import (
    combine_chromosome_stats_tables,
    combine_sim_stats_tables,
)


TABLE_NAMES = [
    "ancestry",
    "kinship",
    "kinship_unrelated",
    "pi_theta_stats",
    "sfs",
    "sfs_2d",
    "ld_decay",
]


parser = argparse.ArgumentParser()
parser.add_argument("--stats-dir", required=True)
parser.add_argument("--num-reps", type=int, required=True)
parser.add_argument("--chroms", nargs="+")
parser.add_argument("--chrom-index", type=int)
parser.add_argument("--chromosomes", action="store_true")
parser.add_argument(
    "--tables",
    nargs="+",
    choices=TABLE_NAMES,
    default=TABLE_NAMES,
)


if __name__ == "__main__":
    args = parser.parse_args()
    if args.chromosomes:
        combine_chromosome_stats_tables(args, args.tables)
    else:
        combine_sim_stats_tables(args, args.tables)

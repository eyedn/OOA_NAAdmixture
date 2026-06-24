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
from python_utils.stats_utils.combine_sim_stats_tables import (
    combine_sim_stats_tables,
)


parser = argparse.ArgumentParser()
parser.add_argument("--stats-dir", required=True)
parser.add_argument("--num-reps", type=int, required=True)
TABLE_NAMES = [
    "ancestry",
    "kinship",
    "pi_theta_stats",
    "sfs",
    "sfs_2d",
    "ld_decay",
]


if __name__ == "__main__":
    combine_sim_stats_tables(parser.parse_args(), TABLE_NAMES)

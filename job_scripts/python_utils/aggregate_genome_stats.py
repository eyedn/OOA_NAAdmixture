###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           aggregate_genome_stats.py
###############################################################################
# aggregate chromosome handoffs into simulation genome statistics.


##### set up ##################################################################
import argparse
from stats_utils.aggregate_genome_stats import (
    aggregate_genome_stats,
)


##### arguments ###############################################################
'''
define command-line arguments for the replicate identifier, input/output
directories, genetic map, chromosomes, and populations.
'''
parser = argparse.ArgumentParser()
parser.add_argument("--rep", type=int, required=True)
parser.add_argument("--stats-dir", required=True)
parser.add_argument("--admixture-dir", required=True)
parser.add_argument("--king-dir", required=True)
parser.add_argument("--genetic-map", required=True)
parser.add_argument("--chroms", nargs="+", required=True)
parser.add_argument("--pops", nargs="+", required=True)


##### main ####################################################################
if __name__ == "__main__":
    args = parser.parse_args()
    aggregate_genome_stats(args)

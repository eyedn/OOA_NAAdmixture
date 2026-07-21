###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           calc_sim_stats.py
###############################################################################

# overview: calculate simulated chromosome statistics from CLI inputs.


##### set up ##################################################################
import argparse
from stats_utils.calc_stats import calc_stats


##### arguments ###############################################################
'''
define command-line arguments for
    - replicate id
    - file paths for ts, plink bed, sample metadata, final ADMIXTURE .fam
    - directory paths for ADMIXTURE, ancestry, KING, and statistics outputs
    - sample size, chrom, genetic map, mutation rate, and list of populations
'''
parser = argparse.ArgumentParser()
parser.add_argument("--rep", type=int, required=True)
parser.add_argument("--tree-tsz-path", required=True)
parser.add_argument("--plink-bed-prefix", required=True)
parser.add_argument("--sample-metadata-path", required=True)
parser.add_argument("--admixture-fam-path", required=True)
parser.add_argument("--admixture-dir", required=True)
parser.add_argument("--global-anc-dir", required=True)
parser.add_argument("--king-dir", required=True)
parser.add_argument("--stats-dir", required=True)
parser.add_argument("--sample-size", type=int, required=True)
parser.add_argument("--chr", required=True)
parser.add_argument("--genetic-map", required=True)
parser.add_argument("--mutation-rate", type=float, required=True)
parser.add_argument("--pops", nargs="+", required=True)


##### main ####################################################################
if __name__ == "__main__":
    args = parser.parse_args()
    calc_stats(args)

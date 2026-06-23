###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           calc_stats.py
###############################################################################


import argparse
from python_utils.stats_utils.calc_stats import calc_stats


parser = argparse.ArgumentParser()
parser.add_argument("--rep", type=int, required=True)
parser.add_argument("--tree-tsz-path", required=True)
parser.add_argument("--plink-bed-prefix", required=True)
parser.add_argument("--admixture-dir", required=True)
parser.add_argument("--king-dir", required=True)
parser.add_argument("--stats-dir", required=True)
parser.add_argument("--sample-size", type=int, required=True)
parser.add_argument("--chr", required=True)
parser.add_argument("--genetic-map", required=True)
parser.add_argument("--mutation-rate", type=float, required=True)
parser.add_argument("--pops", nargs="+", required=True)


if __name__ == "__main__":
    calc_stats(parser.parse_args())

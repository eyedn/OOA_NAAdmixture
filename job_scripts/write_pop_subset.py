###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           write_pop_subset.py
###############################################################################


import argparse
from python_utils.stats_utils.subset_pairs import subset_pairs


# pattern: Imperative Shell


parser = argparse.ArgumentParser()
parser.add_argument("--subset-path", required=True)
parser.add_argument("--pop", required=True)
parser.add_argument("--sample-size", type=int, required=True)
parser.add_argument("--pops", nargs="+", required=True)
args = parser.parse_args()

if __name__ == "__main__":
    subset_pairs(args.subset_path, args.pop, args.pops, args.sample_size)

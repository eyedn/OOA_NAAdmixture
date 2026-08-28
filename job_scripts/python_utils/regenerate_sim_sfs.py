###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           regenerate_sim_sfs.py
###############################################################################

import argparse
from pathlib import Path

import tszip

from shared_utils.write_stats_table import write_stats_table
from sim_utils.calc_stats import _build_1d_sfs_rows


parser = argparse.ArgumentParser()
parser.add_argument("--rep", type=int, required=True)
parser.add_argument("--chrom", required=True)
parser.add_argument("--tree-tsz-path", required=True)
parser.add_argument("--stats-dir", required=True)
parser.add_argument("--sample-size", type=int, required=True)
parser.add_argument("--sfs-size", type=int, required=True)
parser.add_argument("--pops", nargs="+", required=True)


if __name__ == "__main__":
    args = parser.parse_args()
    tree_path = Path(args.tree_tsz_path)
    if not tree_path.is_file() or not str(tree_path).endswith(".ts.tsz"):
        raise FileNotFoundError(f"Missing compressed tree sequence {tree_path}")
    tree_sequence = tszip.decompress(tree_path)
    rows = _build_1d_sfs_rows(
        tree_sequence,
        args.rep,
        args.pops,
        args.sample_size,
        2 * args.sfs_size,
        args.chrom
    )
    output = Path(args.stats_dir) / f"sfs.rep_{args.rep}.chr{args.chrom}"
    write_stats_table(output, rows)

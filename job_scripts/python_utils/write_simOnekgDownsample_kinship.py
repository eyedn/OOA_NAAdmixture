###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           write_simOnekgDownsample_kinship.py
###############################################################################

# overview: combine full-population KING tables for one density-matched task.


##### set up ##################################################################
from pathlib import Path
import argparse

from sim_utils.parse_king_file import parse_king_file
from shared_utils.write_stats_table import write_stats_table


##### arguments ###############################################################
parser = argparse.ArgumentParser()
parser.add_argument("--rep", type=int, required=True)
parser.add_argument("--king-dir", required=True)
parser.add_argument("--stats-dir", required=True)
parser.add_argument("--genetic-map", required=True)
parser.add_argument("--chrom", required=True)
parser.add_argument("--pops", nargs="+", required=True)


##### main ####################################################################
if __name__ == "__main__":
    args = parser.parse_args()
    rows = []
    for pop in args.pops:
        king_path = Path(args.king_dir) / (
            f"{args.genetic_map}_{args.rep}_chr{args.chrom}_{pop}.kin0"
        )
        if not king_path.is_file():
            raise FileNotFoundError(king_path)
        for row in parse_king_file(king_path, args.rep, pop).to_dict(
            "records"
        ):
            rows.append(
                {
                    "rep": row["rep"],
                    "chrom": args.chrom,
                    "pop": row["pop"],
                    "id1": row["sample1"],
                    "id2": row["sample2"],
                    "kinship": row["kinship"]
                }
            )
    write_stats_table(
        Path(args.stats_dir) / f"kinship.rep_{args.rep}.chr{args.chrom}",
        rows,
    )

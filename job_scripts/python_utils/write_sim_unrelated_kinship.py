###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           write_sim_unrelated_kinship.py
###############################################################################

# overview: collect KING unrelated-pair output into the shared kinship schema.


##### set up ##################################################################
from pathlib import Path
import argparse
import csv
import pandas as pd
from sim_utils.parse_king_file import parse_king_file


##### arguments ###############################################################
'''
define command-line arguments for
    - replicate number
    - king and stats output directories
    - genetic map and chromosome
    - list of populations
'''
parser = argparse.ArgumentParser()
parser.add_argument("--rep", type=int, required=True)
parser.add_argument("--king-dir", required=True)
parser.add_argument("--stats-dir", required=True)
parser.add_argument("--genetic-map", required=True)
parser.add_argument("--chr")
parser.add_argument("--pops", nargs="+", required=True)


##### main ####################################################################
if __name__ == "__main__":
    args = parser.parse_args()
    king_dir = Path(args.king_dir)
    stats_dir = Path(args.stats_dir)
    stats_dir.mkdir(parents=True, exist_ok=True)

    # determine whether results are chromosome-level or genome-level
    unit = f"chr{args.chr}" if args.chr is not None else "genome"

    # standardize each population's unrelated-sample KING output.
    rows = []
    for pop in args.pops:
        # create the KING .kin0 path for the unrelated subset
        king_path = king_dir / (
            f"{args.genetic_map}_{args.rep}_{unit}_{pop}_unrelated.kin0"
        )

        # verify that each expected KING file exists
        if not king_path.is_file():
            raise FileNotFoundError(king_path)

        # parse KING output into a standard format using parse_king_file
        king_table = parse_king_file(king_path, args.rep, pop)
        pop_rows = king_table.to_dict("records")
        if args.chr is not None:
            pop_rows = [
                {
                    "rep": row["rep"],
                    "chrom": args.chr,
                    "pop": row["pop"],
                    "sample1": row["sample1"],
                    "sample2": row["sample2"],
                    "kinship": row["kinship"],
                }
                for row in pop_rows
            ]

        # add rows for each pop to a common list of rows
        rows.extend(pop_rows)

    # write the combined unrelated kinship table as a tab-delimited TSV file
    output_suffix = (
        f".rep_{args.rep}.chr{args.chr}"
        if args.chr is not None
        else f".rep_{args.rep}"
    )
    tsv_path = stats_dir / f"kinship_unrelated{output_suffix}.tsv"
    parquet_path = stats_dir / f"kinship_unrelated{output_suffix}.parquet"
    fieldnames = (
        ["rep", "chrom", "pop", "sample1", "sample2", "kinship"]
        if args.chr is not None
        else ["rep", "pop", "sample1", "sample2", "kinship"]
    )
    with open(tsv_path, "w", encoding="utf-8", newline="") as out_file:
        writer = csv.DictWriter(
            out_file,
            fieldnames=fieldnames,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)

    pd.DataFrame(rows, columns=fieldnames).to_parquet(
        parquet_path,
        index=False,
    )

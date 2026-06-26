###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           parse_king_file.py
###############################################################################


import pandas as pd
from .simple_table import SimpleTable


# Parse plink2 KING coefficient output into a table object.
def parse_king_file(king_path, rep, pop):
    if pd is not None:
        king_table = pd.read_csv(king_path, sep=r"\s+")
        sample1_col = "IID1" if "IID1" in king_table.columns else "#IID1"
        sample2_col = "IID2"
        kinship_col = "KINSHIP"
        out_table = pd.DataFrame(
            {
                "rep": rep,
                "pop": pop,
                "sample1": king_table[sample1_col],
                "sample2": king_table[sample2_col],
                "kinship": king_table[kinship_col],
            }
        )
        return out_table

    with open(king_path, "r", encoding="utf-8") as in_file:
        header = in_file.readline().strip().split()
        sample1_col = "IID1" if "IID1" in header else "#IID1"
        sample1_idx = header.index(sample1_col)
        sample2_idx = header.index("IID2")
        kinship_idx = header.index("KINSHIP")
        rows = []
        for line in in_file:
            fields = line.strip().split()
            if not fields:
                continue
            rows.append(
                {
                    "rep": rep,
                    "pop": pop,
                    "sample1": fields[sample1_idx],
                    "sample2": fields[sample2_idx],
                    "kinship": fields[kinship_idx],
                }
            )
    return SimpleTable(rows)

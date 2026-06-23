###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           parse_king_file.py
###############################################################################


# return a dataframe from plink2 KING coefficient calculations
def parse_king_file(king_path, rep, pop):
    import pandas as pd

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

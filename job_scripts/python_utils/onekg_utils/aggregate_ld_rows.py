###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           aggregate_ld_rows.py
###############################################################################
# aggregate chromosome LD-decay summaries at genome scope.


##### set up ##################################################################
from collections import defaultdict
import math


##### main function ###########################################################
'''
pool chromosome LD bins from their summed r2 values and pair counts. Returns
genome-level means without averaging already-averaged chromosome estimates.
'''
def aggregate_ld_rows(rows):
    grouped = defaultdict(lambda: {"sum_r2": 0.0, "n_pairs": 0})
    for row in rows:
        key = (
            int(row["rep"]),
            row["pop"],
            int(row["distance_bin_bp"]),
        )
        grouped[key]["sum_r2"] += float(row["sum_r2"])
        grouped[key]["n_pairs"] += int(row["n_pairs"])
    output = []
    for (rep, pop, distance_bin), values in sorted(grouped.items()):
        num_pairs = values["n_pairs"]
        output.append(
            {
                "rep": rep,
                "pop": pop,
                "distance_bin_bp": distance_bin,
                "mean_r2": (
                    values["sum_r2"] / num_pairs
                    if num_pairs else math.nan
                ),
                "sum_r2": values["sum_r2"],
                "n_pairs": num_pairs,
            }
        )
    return output

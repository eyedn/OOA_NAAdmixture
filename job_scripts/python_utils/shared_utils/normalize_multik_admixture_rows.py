###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           normalize_multik_admixture_rows.py
###############################################################################

# overview: validate unsupervised ancestry rows and preserve component order.

##### set up ##################################################################
import math


MAX_COMPONENTS = 5
Q_SUM_TOLERANCE = 1e-5


##### main function ###########################################################
'''
validate unsupervised ADMIXTURE components and return rows with five neutral
component columns. Components above K are represented by "NA".
'''
def normalize_multik_admixture_rows(rows, k):
    if k < 2 or k > MAX_COMPONENTS:
        raise ValueError(
            f"Ancestry K must be between 2 and {MAX_COMPONENTS}"
        )
    if any(len(row["q_values"]) != k for row in rows):
        raise ValueError(f"Every Q row must contain exactly {k} components")
    for row in rows:
        q_values = row["q_values"]
        if (
            any(not math.isfinite(value) for value in q_values)
            or any(value < 0.0 or value > 1.0 for value in q_values)
            or abs(sum(q_values) - 1.0) > Q_SUM_TOLERANCE
        ):
            raise ValueError("Q rows must contain valid proportions")

    return [
        {
            **row,
            **{
                f"component_{index}_q": (
                    row["q_values"][index - 1] if index <= k else "NA"
                )
                for index in range(1, MAX_COMPONENTS + 1)
            },
        }
        for row in rows
    ]

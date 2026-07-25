###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           orient_multik_admixture_rows.py
###############################################################################

# overview: orient multi-K ADMIXTURE components using AFR/EUR references.

# pattern: Functional Core


##### set up ##################################################################
import math


##### main function ###########################################################
'''
orient unsupervised ADMIXTURE components using African and European reference
populations. Returns rows with stable ancestry-specific columns.
'''
def orient_multik_admixture_rows(rows, afr_pop, eur_pop, k):
    if k < 2:
        raise ValueError("ADMIXTURE K must be at least 2")
    if any(len(row["q_values"]) != k for row in rows):
        raise ValueError(f"Every Q row must contain exactly {k} components")
    for row in rows:
        q_values = row["q_values"]
        if (
            any(not math.isfinite(value) for value in q_values)
            or any(value < 0.0 or value > 1.0 for value in q_values)
            or abs(sum(q_values) - 1.0) > 1e-5
        ):
            raise ValueError("Q rows must contain valid proportions")
    afr_rows = [row for row in rows if row["pop"] == afr_pop]
    eur_rows = [row for row in rows if row["pop"] == eur_pop]
    if not afr_rows or not eur_rows:
        raise ValueError("Both reference populations require Q rows")
    afr_means = [
        sum(row["q_values"][index] for row in afr_rows) / len(afr_rows)
        for index in range(k)
    ]
    afr_index = max(range(k), key=lambda index: afr_means[index])
    eur_means = [
        sum(row["q_values"][index] for row in eur_rows) / len(eur_rows)
        for index in range(k)
    ]
    eur_index = max(
        (index for index in range(k) if index != afr_index),
        key=lambda index: eur_means[index],
    )
    return [
        {
            **row,
            "afr_component": afr_index + 1,
            "eur_component": eur_index + 1,
            "afr_unsupervised_q": row["q_values"][afr_index],
            "eur_unsupervised_q": row["q_values"][eur_index],
            "uncaptured_unsupervised_q": sum(
                value
                for index, value in enumerate(row["q_values"])
                if index not in {afr_index, eur_index}
            ),
        }
        for row in rows
    ]

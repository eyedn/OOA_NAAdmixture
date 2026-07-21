###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           orient_admixture_rows.py
###############################################################################
# orient ADMIXTURE components consistently with African and European labels.


##### main function ###########################################################
'''
orient two unsupervised ADMIXTURE components using African and European
reference populations. Returns rows with stable ancestry-specific columns.
'''
def orient_admixture_rows(rows, afr_pop, eur_pop):
    afr_rows = [row for row in rows if row["pop"] == afr_pop]
    eur_rows = [row for row in rows if row["pop"] == eur_pop]
    if not afr_rows or not eur_rows:
        raise ValueError("Both reference populations require Q rows")
    afr_q1 = sum(row["q1"] for row in afr_rows) / len(afr_rows)
    eur_q1 = sum(row["q1"] for row in eur_rows) / len(eur_rows)
    q1_is_afr = afr_q1 > eur_q1
    return [
        {
            **row,
            "afr_unsupervised_q": row["q1"] if q1_is_afr else row["q2"],
            "eur_unsupervised_q": row["q2"] if q1_is_afr else row["q1"],
        }
        for row in rows
    ]

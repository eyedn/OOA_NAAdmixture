###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           parse_admixture_q.py
###############################################################################


# return a dataframe from the admixture predictions of "admixture --supervised"
def parse_admixture_q(q_path, rep, pops, sample_size):
    import pandas as pd

    rows = []
    q_table = pd.read_csv(q_path, sep=r"\s+", header=None)
    for pop in pops:
        pop_idx = pops.index(pop)
        start = pop_idx * sample_size
        end = (pop_idx + 1) * sample_size
        for sample_offset, row_idx in enumerate(range(start, end), start=1):
            rows.append(
                {
                    "rep": rep,
                    "pop": pop,
                    "sample": f"{pop}_{sample_offset}",
                    "afr_q": q_table.iloc[row_idx, 0],
                    "eur_q": q_table.iloc[row_idx, 1],
                }
            )
    return pd.DataFrame(rows)

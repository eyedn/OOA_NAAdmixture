###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           build_2d_sfs_rows.py
###############################################################################


from itertools import combinations
from .sample_nodes_for_pop import sample_nodes_for_pop


# return the 2d sfs for each pair of populations
def build_2d_sfs_rows(ts, rep, pops, sample_size):
    rows = []
    for pop1, pop2 in combinations(pops, 2):
        sample_sets = [
            sample_nodes_for_pop(ts, pops, sample_size, pop1),
            sample_nodes_for_pop(ts, pops, sample_size, pop2),
        ]
        spectrum = ts.allele_frequency_spectrum(
            sample_sets,
            polarised=True,
            span_normalise=False,
            mode="site",
        )
        for pop1_count in range(spectrum.shape[0]):
            for pop2_count in range(spectrum.shape[1]):
                rows.append(
                    {
                        "rep": rep,
                        "pop1": pop1,
                        "pop2": pop2,
                        "pop1_count": pop1_count,
                        "pop2_count": pop2_count,
                        "count": spectrum[pop1_count, pop2_count],
                    }
                )
    return rows

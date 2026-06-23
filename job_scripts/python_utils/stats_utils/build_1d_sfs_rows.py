###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           build_1d_sfs_rows.py
###############################################################################


from .sample_nodes_for_pop import sample_nodes_for_pop


# return the sfs for each population
def build_1d_sfs_rows(ts, rep, pops, sample_size):
    rows = []
    for pop in pops:
        sample_nodes = sample_nodes_for_pop(ts, pops, sample_size, pop)
        spectrum = ts.allele_frequency_spectrum(
            [sample_nodes],
            polarised=True,
            span_normalise=False,
            mode="site",
        )
        for derived_count, count in enumerate(spectrum):
            rows.append(
                {
                    "rep": rep,
                    "pop": pop,
                    "derived_allele_count": derived_count,
                    "count": count,
                }
            )
    return rows

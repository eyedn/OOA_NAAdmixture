###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           build_pi_theta_rows.py
###############################################################################


from .sample_nodes_for_pop import sample_nodes_for_pop


# pattern: Functional Core


# generate pi and theta results from tree sequence
def build_pi_theta_rows(ts, rep, pops, sample_size, mutation_rate):
    rows = []
    # define watterson's theta as the number of segregating sites of waterson's 
    # constant, normalized to contig length
    for pop in pops:
        sample_nodes = sample_nodes_for_pop(ts, pops, sample_size, pop)
        wattersons_const = sum(
            1.0 / value for value in range(1, len(sample_nodes))
        )
        # returns the number of segregating sites, normalized to contig length
        segregating_sites = ts.segregating_sites(
            [sample_nodes],
            mode="site",
        )[0]

        # append pi and theta values
        pi_value = ts.diversity([sample_nodes], mode="site")[0]
        theta_value = segregating_sites / wattersons_const
        rows.append(
            {
                "rep": rep,
                "pop": pop,
                "stat": "pi",
                "value": pi_value,
                "ne_value": pi_value / (4 * mutation_rate),
                "mutation_rate": mutation_rate,
            }
        )
        rows.append(
            {
                "rep": rep,
                "pop": pop,
                "stat": "theta",
                "value": theta_value,
                "ne_value": theta_value / (4 * mutation_rate),
                "mutation_rate": mutation_rate,
            }
        )
    return rows

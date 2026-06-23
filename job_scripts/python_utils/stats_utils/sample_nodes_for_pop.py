###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           sample_nodes_for_pop.py
###############################################################################


# return all nodes in the tree sequnce for a specific pop
def sample_nodes_for_pop(ts, pops, sample_size, pop):
    pop_idx = pops.index(pop)
    start = pop_idx * sample_size
    end = (pop_idx + 1) * sample_size
    nodes = []
    for ind_id in range(start, end):
        nodes.extend(int(node) for node in ts.individual(ind_id).nodes)
    return nodes

###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           get_local_ancestry_table.py
###############################################################################


import tspop


# Get local ancestry via the census before the first admixture pulse.
def get_local_ancestry_table(ts, census_time):
    return tspop.get_pop_ancestry(
        ts,
        census_time=census_time,
    ).ancestry_table.copy()

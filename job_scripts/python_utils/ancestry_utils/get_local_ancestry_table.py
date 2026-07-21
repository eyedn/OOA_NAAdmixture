###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           get_local_ancestry_table.py
###############################################################################

# overview: extract per-sample local ancestry intervals from tree sequence.


##### set up ##################################################################
import tspop


##### main function ###########################################################
'''
extract local ancestry intervals through the census placed before the first
admixture pulse. Returns a copy of the tspop ancestry table.
'''
def get_local_ancestry_table(ts, census_time):
    return tspop.get_pop_ancestry(
        ts,
        census_time=census_time,
    ).ancestry_table.copy()

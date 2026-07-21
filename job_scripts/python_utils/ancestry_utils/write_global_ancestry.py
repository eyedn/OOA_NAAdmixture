###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           write_global_ancestry.py
###############################################################################
# overview: write global ancestry tables for downstream statistics.


##### main function ###########################################################
'''
write derived global ancestry proportions to a tab-separated table with column
names and without a pandas index.
'''
def write_global_ancestry(global_ancestry_table, out_tsv):
    global_ancestry_table.to_csv(
        out_tsv,
        sep="\t",
        header=True,
        index=False,
    )

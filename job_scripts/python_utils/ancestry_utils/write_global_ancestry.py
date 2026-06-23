###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           write_global_ancestry.py
###############################################################################


# output tsv of derived global ancestry
def write_global_ancestry(global_ancestry_table, out_tsv):
    global_ancestry_table.to_csv(
        out_tsv,
        sep="\t",
        header=True,
        index=False,
    )

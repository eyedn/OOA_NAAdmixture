###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           aggregate_folded_sfs_rows.py
###############################################################################
# aggregate folded one-dimensional SFS rows across chromosome inputs.


##### set up ##################################################################
from collections import defaultdict


##### main function ###########################################################
'''
sum folded one-dimensional SFS bins across chromosome tables. Returns one row
per replicate, population, and minor-allele count.
'''
def aggregate_folded_sfs_rows(rows):
    grouped = defaultdict(float)
    for row in rows:
        key = (
            int(row["rep"]),
            row["pop"],
            int(row["minor_allele_count"]),
        )
        grouped[key] += float(row["count"])
    return [
        {
            "rep": rep,
            "pop": pop,
            "minor_allele_count": minor_count,
            "count": count,
        }
        for (rep, pop, minor_count), count in sorted(grouped.items())
    ]

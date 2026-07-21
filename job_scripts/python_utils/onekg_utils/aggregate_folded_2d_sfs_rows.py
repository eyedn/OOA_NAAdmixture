###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           aggregate_folded_2d_sfs_rows.py
###############################################################################
# aggregate folded two-dimensional SFS rows across chromosome inputs.


##### set up ##################################################################
from collections import defaultdict


##### main function ###########################################################
'''
sum folded pairwise SFS bins across chromosome tables. Returns one row per
replicate, population pair, and pair of minor-allele counts.
'''
def aggregate_folded_2d_sfs_rows(rows):
    grouped = defaultdict(float)
    for row in rows:
        key = (
            int(row["rep"]),
            row["pop1"],
            row["pop2"],
            int(row["pop1_minor_allele_count"]),
            int(row["pop2_minor_allele_count"]),
        )
        grouped[key] += float(row["count"])
    return [
        {
            "rep": rep,
            "pop1": pop1,
            "pop2": pop2,
            "pop1_minor_allele_count": pop1_count,
            "pop2_minor_allele_count": pop2_count,
            "count": count,
        }
        for (
            rep,
            pop1,
            pop2,
            pop1_count,
            pop2_count,
        ), count in sorted(grouped.items())
    ]

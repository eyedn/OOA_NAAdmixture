###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           build_folded_2d_sfs_rows.py
###############################################################################
# build folded two-dimensional SFS rows for population comparisons.


##### set up ##################################################################
from collections import Counter


##### main function ###########################################################
'''
build folded pairwise SFS rows at positions shared by both populations. Returns
minor-allele-count bins for the two populations at their shared sites.
'''
def build_folded_2d_sfs_rows(
    rep,
    chrom,
    pop1,
    pop2,
    pop1_counts,
    pop2_counts,
):
    shared_positions = sorted(set(pop1_counts) & set(pop2_counts))
    bins = Counter(
        (
            min(pop1_counts[position]),
            min(pop2_counts[position]),
        )
        for position in shared_positions
    )
    return [
        {
            "rep": rep,
            "chrom": chrom,
            "pop1": pop1,
            "pop2": pop2,
            "pop1_minor_allele_count": pop1_count,
            "pop2_minor_allele_count": pop2_count,
            "count": count,
        }
        for (pop1_count, pop2_count), count in sorted(bins.items())
    ]

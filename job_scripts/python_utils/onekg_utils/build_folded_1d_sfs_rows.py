###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           build_folded_1d_sfs_rows.py
###############################################################################
# build folded one-dimensional SFS rows without ancestral-state inference.


##### set up ##################################################################
from collections import Counter


##### main function ###########################################################
'''
build a folded one-dimensional SFS from reference and alternate allele counts.
Returns minor-allele-count bins without inferring ancestral state.
'''
def build_folded_1d_sfs_rows(rep, chrom, pop, allele_counts):
    bins = Counter(min(ref_count, alt_count)
                   for ref_count, alt_count in allele_counts)
    return [
        {
            "rep": rep,
            "chrom": chrom,
            "pop": pop,
            "minor_allele_count": minor_count,
            "count": count,
        }
        for minor_count, count in sorted(bins.items())
    ]

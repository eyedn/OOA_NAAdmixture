###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           build_sample_ids.py
###############################################################################

# overview: construct stable sample identifiers for ancestry outputs.


##### main function ###########################################################
'''
build tree-sequence and VCF identifiers for one population sample. Returns the
identifier in full tree-sequence order and the population-relative VCF label.
'''
def build_sample_ids(pop, sample_ind, pop_start_ind):
    sample_id = f"{pop}_{sample_ind + 1}"
    vcf_sample_id = f"{pop}_{sample_ind - pop_start_ind + 1}"

    return sample_id, vcf_sample_id

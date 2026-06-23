###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           build_adx_sample_ids.py
###############################################################################

# pattern: Functional Core


# build both ADX labels for a tree-sequence individual index
def build_adx_sample_ids(sample_ind, adx_start_ind):
    sample_id = f"ADX_{sample_ind + 1}"
    vcf_sample_id = f"ADX_{sample_ind - adx_start_ind + 1}"

    return sample_id, vcf_sample_id

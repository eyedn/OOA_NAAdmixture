###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           build_sample_ids.py
###############################################################################


# build tree-sequence and VCF labels for a population-specific sample
def build_sample_ids(pop, sample_ind, pop_start_ind):
    sample_id = f"{pop}_{sample_ind + 1}"
    vcf_sample_id = f"{pop}_{sample_ind - pop_start_ind + 1}"

    return sample_id, vcf_sample_id

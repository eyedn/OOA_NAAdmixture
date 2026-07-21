###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           read_sample_metadata.py
###############################################################################
# read simulation metadata used to build supervised ADMIXTURE labels.


##### set up ##################################################################
import csv


##### main function ###########################################################
'''
read sample metadata rows for supervised ADMIXTURE labels. Returns all metadata
rows as a list for downstream supervised ADMIXTURE processing.
'''
def read_sample_metadata(sample_metadata_path):
    with open(sample_metadata_path, "r", encoding="utf-8") as in_file:
        return list(csv.DictReader(in_file, delimiter="\t"))

###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           read_tsv_rows.py
###############################################################################
# read canonical tab-separated statistics rows.


##### set up ##################################################################
import csv


##### main function ###########################################################
'''
read one tab-separated file as a list of dictionaries while preserving the
column names and row order.
'''
def read_tsv_rows(path):
    with open(path, "r", encoding="utf-8", newline="") as in_file:
        return list(csv.DictReader(in_file, delimiter="\t"))

###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           read_plink_ld_rows.py
###############################################################################

# overview: read PLINK LD output rows for distance-bin aggregation.


##### main function ###########################################################
'''
stream dictionary rows from a whitespace-delimited PLINK LD report. Header
names are normalized by removing PLINK's optional leading hash character.
'''
def read_plink_ld_rows(path):
    with open(path, "r", encoding="utf-8") as in_file:
        header_line = next((line for line in in_file if line.strip()), None)
        if header_line is None:
            raise ValueError(f"Empty PLINK LD output {path}")
        headers = [header.lstrip("#") for header in header_line.split()]
        for line in in_file:
            if line.strip():
                yield dict(zip(headers, line.split()))

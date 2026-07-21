###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           parse_gt.py
###############################################################################
# parse VCF genotype fields into called allele counts.


##### set up ##################################################################
import re


##### main function ###########################################################
'''
convert one diploid GT field to reference and alternate allele counts. Returns
None for missing genotypes and rejects non-biallelic or non-diploid values.
'''
def parse_gt(gt):
    alleles = re.split(r"[|/]", gt)
    if len(alleles) != 2:
        raise ValueError(f"Expected diploid GT, found {gt}")
    if any(allele == "." for allele in alleles):
        return None
    if any(allele not in {"0", "1"} for allele in alleles):
        raise ValueError(f"Invalid biallelic allele index in GT {gt}")
    alt_count = sum(allele == "1" for allele in alleles)
    return 2 - alt_count, alt_count

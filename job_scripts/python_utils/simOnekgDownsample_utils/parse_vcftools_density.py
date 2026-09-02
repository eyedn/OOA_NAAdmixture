###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           parse_vcftools_density.py
###############################################################################

# overview: normalize vcftools SNP-density rows for density selection.


##### main ####################################################################
'''
parse vcftools --SNPdensity output from an iterable of text lines.
'''
def parse_vcftools_density(lines):
    nonempty = [line.rstrip("\n").split("\t") for line in lines if line.strip()]
    if not nonempty:
        raise ValueError("SNP-density table is empty")
    headers = [value.upper() for value in nonempty[0]]
    required = ("CHROM", "BIN_START", "SNP_COUNT")
    if any(name not in headers for name in required):
        raise ValueError("SNP-density table lacks required vcftools columns")
    indexes = {name: headers.index(name) for name in required}
    rows = []
    for fields in nonempty[1:]:
        if max(indexes.values()) >= len(fields):
            raise ValueError("Malformed SNP-density row")
        rows.append(
            {
                "chrom": fields[indexes["CHROM"]],
                "bin_start": int(fields[indexes["BIN_START"]]),
                "snp_count": int(fields[indexes["SNP_COUNT"]])
            }
        )
    return rows

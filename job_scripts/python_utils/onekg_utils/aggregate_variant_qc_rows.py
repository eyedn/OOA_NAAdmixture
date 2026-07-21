###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           aggregate_variant_qc_rows.py
###############################################################################
# aggregate chromosome variant-quality counts at genome scope.


##### main function ###########################################################
'''
sum chromosome variant-QC counts while retaining population sample counts.
Returns one genome-level row for the empirical replicate.
'''
def aggregate_variant_qc_rows(rows):
    output = {"rep": 0, "chrom": "genome", "pop": "ALL"}
    for row in rows:
        for key, value in row.items():
            if key in {"rep", "chrom", "pop"} or value in {None, ""}:
                continue
            number = int(float(value))
            if key.startswith("retained_") and key.endswith("_samples"):
                output[key] = max(output.get(key, 0), number)
            else:
                output[key] = output.get(key, 0) + number
    return [output]

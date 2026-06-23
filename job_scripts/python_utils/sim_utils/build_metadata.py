###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           build_metadata.py
###############################################################################


# generate metadata for "admixture --supervised"
def build_metadata(pops, sample_size):
    rows = []
    original_order = 1
    for pop in pops:
        supervised_label = "-" if pop == "ADX" else pop
        for sample_idx in range(1, sample_size + 1):
            rows.append(
                {
                    "fid": "0",
                    "iid": f"{pop}_{sample_idx}",
                    "pop": pop,
                    "supervised_label": supervised_label,
                    "original_order": original_order,
                }
            )
            original_order += 1
    return rows

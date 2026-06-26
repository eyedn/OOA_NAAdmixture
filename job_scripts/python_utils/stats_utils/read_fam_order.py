###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           read_fam_order.py
###############################################################################


# read FID and IID values in the order retained by a PLINK FAM file
def read_fam_order(fam_path):
    rows = []
    with open(fam_path, "r", encoding="utf-8") as in_file:
        for line in in_file:
            fields = line.strip().split()
            if not fields:
                continue
            rows.append({"fid": fields[0], "iid": fields[1]})
    return rows


# read ADMIXTURE Q rows using the filtered FAM sample order
def read_q_rows(rep, q_path, fam_path):
    fam_table = read_fam_order(fam_path)
    q_rows = []
    with open(q_path, "r", encoding="utf-8") as in_file:
        for line in in_file:
            fields = line.strip().split()
            if not fields:
                continue
            q_rows.append((float(fields[0]), float(fields[1])))

    if len(q_rows) != len(fam_table):
        raise ValueError(
            "ADMIXTURE Q row count does not match final FAM row count"
        )

    rows = []
    for fam_row, q_row in zip(fam_table, q_rows):
        pop = fam_row["iid"].split("_", 1)[0]
        rows.append(
            {
                "rep": rep,
                "pop": pop,
                "vcf_sample_id": fam_row["iid"],
                "afr_q": q_row[0],
                "eur_q": q_row[1],
            }
        )

    return rows

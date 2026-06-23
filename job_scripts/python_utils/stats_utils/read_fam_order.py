###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           read_fam_order.py
###############################################################################


# pattern: Imperative Shell


def read_fam_order(fam_path):
    rows = []
    with open(fam_path, "r", encoding="utf-8") as in_file:
        for line in in_file:
            fields = line.strip().split()
            if not fields:
                continue
            rows.append({"fid": fields[0], "iid": fields[1]})
    return rows

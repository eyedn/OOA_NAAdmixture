###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           parse_admixture_q.py
###############################################################################


# pattern: Imperative Shell


# return a dataframe from the admixture predictions of "admixture --supervised"
def parse_admixture_q(q_path, rep, sample_metadata_path, fam_path):
    from .read_fam_order import read_fam_order
    from .read_sample_metadata import read_sample_metadata
    from .simple_table import SimpleTable

    fam_table = read_fam_order(fam_path)
    metadata_by_sample = {
        (row["fid"], row["iid"]): row
        for row in read_sample_metadata(sample_metadata_path)
    }
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
        sample_key = (fam_row["fid"], fam_row["iid"])
        if sample_key not in metadata_by_sample:
            raise ValueError(
                f"Missing sample metadata for: {fam_row['iid']}"
            )
        metadata_row = metadata_by_sample[sample_key]
        rows.append(
            {
                "rep": rep,
                "pop": metadata_row["pop"],
                "sample": fam_row["iid"],
                "afr_q": q_row[0],
                "eur_q": q_row[1],
            }
        )

    try:
        import pandas as pd
    except (ImportError, ValueError):
        return SimpleTable(rows)
    return pd.DataFrame(rows)

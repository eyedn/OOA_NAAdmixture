###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           build_ancestry_table.py
###############################################################################

# pattern: Imperative Shell


import csv
from pathlib import Path


ANCESTRY_COLUMNS = [
    "rep",
    "pop",
    "sample_id",
    "vcf_sample_id",
    "afr_tspop",
    "eur_tspop",
    "afr_q",
    "eur_q",
]


def read_q_rows(rep, q_path, fam_path):
    from .read_fam_order import read_fam_order

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


# build joined tspop and ADMIXTURE ancestry rows for one replicate
def build_ancestry_rows(
    rep,
    pops,
    genetic_map,
    global_anc_dir,
    q_path,
    fam_path,
):
    global_rows = []
    for pop in pops:
        global_path = Path(global_anc_dir) / f"{genetic_map}_{rep}_{pop}.tsv"
        with open(global_path, "r", encoding="utf-8") as in_file:
            reader = csv.DictReader(in_file, delimiter="\t")
            for row in reader:
                global_rows.append(
                    {
                        "rep": rep,
                        "pop": pop,
                        "sample_id": row["sample_id"],
                        "vcf_sample_id": row["vcf_sample_id"],
                        "afr_tspop": float(row["AFR_prop"]),
                        "eur_tspop": float(row["EUR_prop"]),
                    }
                )

    q_by_sample = {
        (row["pop"], row["vcf_sample_id"]): row
        for row in read_q_rows(rep, q_path, fam_path)
    }
    ancestry_rows = []
    for row in global_rows:
        q_row = q_by_sample.get((row["pop"], row["vcf_sample_id"]), {})
        ancestry_rows.append(
            {
                **row,
                "afr_q": q_row.get("afr_q"),
                "eur_q": q_row.get("eur_q"),
            }
        )

    return ancestry_rows


# build joined tspop and ADMIXTURE ancestry table for one replicate
def build_ancestry_table(*args, **kwargs):
    from .simple_table import SimpleTable

    rows = build_ancestry_rows(*args, **kwargs)
    try:
        import pandas as pd
    except (ImportError, ValueError):
        return SimpleTable(rows)
    return pd.DataFrame(rows, columns=ANCESTRY_COLUMNS)

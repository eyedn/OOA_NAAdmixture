###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           build_onekg_ancestry.py
###############################################################################

# overview: build empirical ancestry tables from sample labels and ADMIXTURE
# outputs.



##### set up ##################################################################
from pathlib import Path
import argparse
from onekg_utils.orient_admixture_rows import orient_admixture_rows
from onekg_utils.normalize_multik_admixture_rows import (
    normalize_multik_admixture_rows,
)
from onekg_utils.read_onekg_sample_pops import read_onekg_sample_pops
from onekg_utils.write_stats_table import write_stats_table


##### arguments ###############################################################
'''
define command-line arguments for empirical sample metadata, supervised and
unsupervised Q files, population labels, and the ancestry output directory.
'''
parser = argparse.ArgumentParser()
parser.add_argument("--unrels-path", required=True)
parser.add_argument("--source-fam-path", required=True)
parser.add_argument("--admixture-fam-path", required=True)
parser.add_argument("--supervised-q-path", required=True)
parser.add_argument(
    "--unsupervised-q-path",
    action="append",
    required=True,
    metavar="K=PATH",
)
parser.add_argument("--stats-dir", required=True)
parser.add_argument("--chrom")
parser.add_argument("--afr-pop", required=True)
parser.add_argument("--eur-pop", required=True)
parser.add_argument("--admixed-pop", required=True)


##### main ####################################################################
if __name__ == "__main__":
    args = parser.parse_args()
    pops = [args.afr_pop, args.eur_pop, args.admixed_pop]
    sample_pops = read_onekg_sample_pops(
        args.unrels_path,
        args.source_fam_path,
        pops,
    )
    with open(args.admixture_fam_path, "r", encoding="utf-8") as in_file:
        samples = [line.split()[1] for line in in_file if line.strip()]

    # read and orient supervised K=2 estimates in final FAM order.
    with open(args.supervised_q_path, "r", encoding="utf-8") as in_file:
        supervised_values = [
            line.split() for line in in_file if line.strip()
        ]
    if len(supervised_values) != len(samples):
        raise ValueError(
            f"Q/FAM row count mismatch for {args.supervised_q_path}"
        )
    supervised_rows = [
        {
            "sample": sample,
            "pop": sample_pops[sample],
            "q1": float(q_values[0]),
            "q2": float(q_values[1]),
        }
        for sample, q_values in zip(samples, supervised_values)
    ]
    supervised = orient_admixture_rows(
        supervised_rows,
        args.afr_pop,
        args.eur_pop,
    )

    # map each unsupervised K to its explicitly preserved Q-file path.
    unsupervised_q_paths = {}
    for q_spec in args.unsupervised_q_path:
        k_text, separator, q_path = q_spec.partition("=")
        if not separator or not k_text.isdigit() or not q_path:
            raise ValueError(
                "Unsupervised Q paths must use the format K=PATH"
            )
        k = int(k_text)
        if k in unsupervised_q_paths:
            raise ValueError(f"Duplicate unsupervised ADMIXTURE K={k}")
        unsupervised_q_paths[k] = q_path
    if 2 not in unsupervised_q_paths:
        raise ValueError("Unsupervised ADMIXTURE K=2 is required")

    # preserve every unsupervised matrix in K-then-FAM row order.
    multik_rows = []
    for k, q_path in sorted(unsupervised_q_paths.items()):
        with open(q_path, "r", encoding="utf-8") as in_file:
            values = [line.split() for line in in_file if line.strip()]
        if len(values) != len(samples):
            raise ValueError(f"Q/FAM row count mismatch for {q_path}")
        raw_rows = [
            {
                "sample": sample,
                "pop": sample_pops[sample],
                "q_values": [float(value) for value in q_values],
            }
            for sample, q_values in zip(samples, values)
        ]
        normalized = normalize_multik_admixture_rows(raw_rows, k)
        for normalized_row in normalized:
            row = {"rep": 0}
            if args.chrom is not None:
                row["chrom"] = args.chrom
            row.update(
                {
                    "pop": normalized_row["pop"],
                    "sample_id": "NA",
                    "vcf_sample_id": normalized_row["sample"],
                    "k": k,
                    "component_1_q": normalized_row["component_1_q"],
                    "component_2_q": normalized_row["component_2_q"],
                    "component_3_q": normalized_row["component_3_q"],
                    "component_4_q": normalized_row["component_4_q"],
                    "component_5_q": normalized_row["component_5_q"],
                    "span": "NA",
                }
            )
            multik_rows.append(row)

    # format supervised estimates in the canonical ancestry table schema.
    rows = []
    for supervised_row in supervised:
        row = {"rep": 0}
        if args.chrom is not None:
            row["chrom"] = args.chrom
        row.update(
            {
                "pop": supervised_row["pop"],
                "sample_id": "NA",
                "vcf_sample_id": supervised_row["sample"],
                "afr_tspop": "NA",
                "eur_tspop": "NA",
                "afr_q": supervised_row["afr_unsupervised_q"],
                "eur_q": supervised_row["eur_unsupervised_q"],
                "span": "NA",
            }
        )
        rows.append(row)
    suffix = f".chr{args.chrom}" if args.chrom is not None else ""
    output = Path(args.stats_dir) / f"ancestry.rep_0{suffix}"
    write_stats_table(output, rows)
    multik_output = (
        Path(args.stats_dir) / f"ancestry_multik.rep_0{suffix}"
    )
    write_stats_table(multik_output, multik_rows)

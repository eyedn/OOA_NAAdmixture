###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           build_onekg_ancestry.py
###############################################################################
# build empirical ancestry tables from sample labels and ADMIXTURE outputs.


##### set up ##################################################################
from pathlib import Path
import argparse
from onekg_utils.orient_admixture_rows import orient_admixture_rows
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
parser.add_argument("--unsupervised-q-path", required=True)
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
    # read and orient supervised and unsupervised Q matrices in FAM order.
    q_sets = []
    for q_path in (args.supervised_q_path, args.unsupervised_q_path):
        with open(q_path, "r", encoding="utf-8") as in_file:
            values = [line.split() for line in in_file if line.strip()]
        if len(values) != len(samples):
            raise ValueError(f"Q/FAM row count mismatch for {q_path}")
        raw_rows = [
            {
                "sample": sample,
                "pop": sample_pops[sample],
                "q1": float(q_values[0]),
                "q2": float(q_values[1]),
            }
            for sample, q_values in zip(samples, values)
        ]
        q_sets.append(
            orient_admixture_rows(raw_rows, args.afr_pop, args.eur_pop)
        )
    supervised, unsupervised = q_sets
    # join both ADMIXTURE results into the canonical ancestry table schema.
    rows = []
    for supervised_row, unsupervised_row in zip(supervised, unsupervised):
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
                "afr_unsupervised_q": (
                    unsupervised_row["afr_unsupervised_q"]
                ),
                "eur_unsupervised_q": (
                    unsupervised_row["eur_unsupervised_q"]
                ),
                "span": "NA",
            }
        )
        rows.append(row)
    suffix = f".chr{args.chrom}" if args.chrom is not None else ""
    output = Path(args.stats_dir) / f"ancestry.rep_0{suffix}"
    write_stats_table(output, rows)

###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           build_onekg_ancestry.py
###############################################################################

# overview: build empirical ADMIXTURE and fastStructure ancestry tables.



##### set up ##################################################################
from pathlib import Path
import argparse
from onekg_utils.build_multik_ancestry_rows import (
    build_multik_ancestry_rows,
)
from onekg_utils.orient_admixture_rows import orient_admixture_rows
from onekg_utils.parse_faststructure_choose_k import (
    parse_faststructure_choose_k,
)
from onekg_utils.parse_k_path_specs import parse_k_path_specs
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
    "--admixture-q-path",
    action="append",
    required=True,
    metavar="K=PATH",
)
parser.add_argument(
    "--faststructure-q-path",
    action="append",
    required=True,
    metavar="K=PATH",
)
parser.add_argument("--faststructure-choose-k-path", required=True)
parser.add_argument("--faststructure-prior", required=True)
parser.add_argument("--faststructure-seed", required=True, type=int)
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

    # read both unsupervised tools against one authoritative final FAM order.
    admixture_q_paths = parse_k_path_specs(
        args.admixture_q_path,
        "ADMIXTURE",
    )
    faststructure_q_paths = parse_k_path_specs(
        args.faststructure_q_path,
        "fastStructure",
    )
    if sorted(admixture_q_paths) != sorted(faststructure_q_paths):
        raise ValueError(
            "ADMIXTURE and fastStructure must use the same K values"
        )
    values_by_tool = {}
    for tool_name, paths in (
        ("ADMIXTURE", admixture_q_paths),
        ("fastStructure", faststructure_q_paths),
    ):
        values_by_k = {}
        for k, q_path in sorted(paths.items()):
            with open(q_path, "r", encoding="utf-8") as in_file:
                values_by_k[k] = [
                    [float(value) for value in line.split()]
                    for line in in_file
                    if line.strip()
                ]
        values_by_tool[tool_name] = values_by_k
    admixture_multik_rows = build_multik_ancestry_rows(
        samples,
        sample_pops,
        values_by_tool["ADMIXTURE"],
        args.chrom,
    )
    faststructure_multik_rows = build_multik_ancestry_rows(
        samples,
        sample_pops,
        values_by_tool["fastStructure"],
        args.chrom,
    )

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
    stats_dir = Path(args.stats_dir)
    write_stats_table(
        stats_dir / f"ancestry_ADMIXTURE_super.rep_0{suffix}",
        rows,
    )
    write_stats_table(
        stats_dir / f"ancestry_ADMIXTURE_multik.rep_0{suffix}",
        admixture_multik_rows,
    )
    write_stats_table(
        stats_dir / f"ancestry_fastStructure_multik.rep_0{suffix}",
        faststructure_multik_rows,
    )
    with open(
        args.faststructure_choose_k_path,
        "r",
        encoding="utf-8",
    ) as in_file:
        choose_k_report = in_file.read()
    choose_k_row = parse_faststructure_choose_k(
        choose_k_report,
        args.faststructure_prior,
        args.faststructure_seed,
        sorted(faststructure_q_paths),
        args.chrom,
    )
    write_stats_table(
        stats_dir / f"fastStructure_chooseK.rep_0{suffix}",
        [choose_k_row],
    )

    # remove superseded empirical filenames after all replacement tables exist.
    for old_family in ("ancestry", "ancestry_multik"):
        old_prefix = stats_dir / f"{old_family}.rep_0{suffix}"
        for extension in ("tsv", "parquet"):
            Path(f"{old_prefix}.{extension}").unlink(missing_ok=True)

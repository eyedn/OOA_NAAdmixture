###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           build_sim_inference_ancestry.py
###############################################################################

# overview: build simulated unsupervised ancestry and choose-K tables.


##### set up ##################################################################
from pathlib import Path
import argparse
import csv
from shared_utils.build_multik_ancestry_rows import build_multik_ancestry_rows
from shared_utils.parse_faststructure_choose_k import parse_faststructure_choose_k
from shared_utils.parse_k_path_specs import parse_k_path_specs
from shared_utils.write_stats_table import write_stats_table


##### arguments ###############################################################
parser = argparse.ArgumentParser()
parser.add_argument("--rep", required=True, type=int)
parser.add_argument("--sample-metadata-path", required=True)
parser.add_argument("--admixture-fam-path", required=True)
parser.add_argument(
    "--admixture-q-path",
    action="append",
    required=True,
    metavar="K=PATH"
)
parser.add_argument(
    "--faststructure-q-path",
    action="append",
    required=True,
    metavar="K=PATH"
)
parser.add_argument("--faststructure-choose-k-path", required=True)
parser.add_argument("--faststructure-prior", required=True)
parser.add_argument("--faststructure-seed", required=True, type=int)
parser.add_argument("--stats-dir", required=True)
parser.add_argument("--chrom")


##### main ####################################################################
if __name__ == "__main__":
    args = parser.parse_args()
    with open(
        args.sample_metadata_path,
        "r",
        encoding="utf-8",
        newline=""
    ) as in_file:
        metadata_rows = list(csv.DictReader(in_file, delimiter="\t"))
    metadata_by_sample = {row["iid"]: row for row in metadata_rows}

    with open(args.admixture_fam_path, "r", encoding="utf-8") as in_file:
        samples = [line.split()[1] for line in in_file if line.strip()]
    missing_samples = [
        sample for sample in samples if sample not in metadata_by_sample
    ]
    if missing_samples:
        raise ValueError(
            "Missing simulation metadata for FAM samples: "
            + ", ".join(missing_samples)
        )
    sample_pops = {
        sample: metadata_by_sample[sample]["pop"] for sample in samples
    }
    sample_ids = {
        sample: (
            f"{metadata_by_sample[sample]['pop']}_"
            f"{metadata_by_sample[sample]['original_order']}"
        )
        for sample in samples
    }

    admixture_paths = parse_k_path_specs(
        args.admixture_q_path,
        "ADMIXTURE"
    )
    faststructure_paths = parse_k_path_specs(
        args.faststructure_q_path,
        "fastStructure"
    )
    if sorted(admixture_paths) != sorted(faststructure_paths):
        raise ValueError(
            "ADMIXTURE and fastStructure must use the same K values"
        )

    ancestry_by_tool = {}
    for table_name, paths in (
        ("ancestry_ADMIXTURE_multik", admixture_paths),
        ("ancestry_fastStructure_multik", faststructure_paths)
    ):
        values_by_k = {}
        for k, q_path in sorted(paths.items()):
            with open(q_path, "r", encoding="utf-8") as in_file:
                values_by_k[k] = [
                    [float(value) for value in line.split()]
                    for line in in_file
                    if line.strip()
                ]
        ancestry_by_tool[table_name] = values_by_k
    suffix = f".chr{args.chrom}" if args.chrom is not None else ""
    stats_dir = Path(args.stats_dir)
    for table_name, values_by_k in ancestry_by_tool.items():
        rows = build_multik_ancestry_rows(
            samples,
            sample_pops,
            values_by_k,
            chrom=args.chrom,
            rep=args.rep,
            sample_ids=sample_ids
        )
        write_stats_table(
            stats_dir / f"{table_name}.rep_{args.rep}{suffix}",
            rows,
        )

    with open(
        args.faststructure_choose_k_path,
        "r",
        encoding="utf-8"
    ) as in_file:
        choose_k_report = in_file.read()
    choose_k_row = parse_faststructure_choose_k(
        choose_k_report,
        args.faststructure_prior,
        args.faststructure_seed,
        sorted(faststructure_paths),
        chrom=args.chrom,
        rep=args.rep
    )
    write_stats_table(
        stats_dir / f"fastStructure_chooseK.rep_{args.rep}{suffix}",
        [choose_k_row]
    )

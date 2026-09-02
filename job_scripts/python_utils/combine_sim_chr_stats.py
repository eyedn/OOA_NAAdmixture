###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           combine_sim_chr_stats.py
###############################################################################

# overview: combine simulated chromosome statistics across replicate outputs.


##### set up ##################################################################
from pathlib import Path
import argparse
import csv
import subprocess
import sys

from combine_sim_chr_stats_utils.validate_parquet_table import (
    validate_parquet_table,
)

'''
define supported chrom-level table types:
    - ancestry
    - kinship
    - unrelated kinship
    - pi/theta
    - one-dimensional SFS
    - LD decay
'''
TABLE_NAMES = [
    "ancestry",
    "ancestry_ADMIXTURE_super",
    "ancestry_ADMIXTURE_multik",
    "ancestry_fastStructure_multik",
    "fastStructure_chooseK",
    "kinship",
    "kinship_unrelated",
    "pi_theta_stats",
    "sfs",
    "ld_decay"
]
SIM_ONEKG_DOWNSAMPLE_TABLE_NAMES = [
    "ancestry_ADMIXTURE_super",
    "ancestry_ADMIXTURE_multik",
    "ancestry_fastStructure_multik",
    "fastStructure_chooseK",
    "kinship",
    "kinship_unrelated",
    "pi_theta_stats_intergenic",
    "pi_theta_stats_full_callable_chrom",
    "sfs",
    "ld_decay",
    "variant_qc",
    "snp_density"
]
TABLE_CHOICES = sorted(set(TABLE_NAMES + SIM_ONEKG_DOWNSAMPLE_TABLE_NAMES))
NONEMPTY_TABLE_NAMES = {
    "ancestry_ADMIXTURE_super",
    "ancestry_ADMIXTURE_multik",
    "ancestry_fastStructure_multik",
    "fastStructure_chooseK",
    "pi_theta_stats_intergenic",
    "pi_theta_stats_full_callable_chrom",
    "sfs",
    "variant_qc",
    "snp_density"
}


##### internal functions #####################################################
'''
internal: validate and concatenate required replicate TSV/Parquet table pairs.
'''
def _combine_table_paths(
    paths,
    expected_reps=None,
    expected_chrom=None,
    require_nonempty=False,
    return_schema=False
):
    paths = [Path(path) for path in paths]
    companion_paths = [path.with_suffix(".parquet") for path in paths]
    missing_paths = [
        path for path in paths + companion_paths if not path.exists()
    ]
    if missing_paths:
        raise FileNotFoundError(", ".join(str(path) for path in missing_paths))
    if expected_reps is not None and len(expected_reps) != len(paths):
        raise ValueError("Expected replicate count does not match input paths")
    rows = []
    expected_schema = None
    for path_index, path in enumerate(paths):
        with open(path, "r", encoding="utf-8", newline="") as in_file:
            reader = csv.DictReader(in_file, delimiter="\t")
            schema = reader.fieldnames
            if not schema:
                raise ValueError(f"Missing TSV schema in {path}")
            if expected_schema is None:
                expected_schema = schema
            elif schema != expected_schema:
                raise ValueError(f"tsv schema mismatch in {path}")
            path_rows = list(reader)
        malformed_rows = [
            row for row in path_rows
            if None in row or any(value is None for value in row.values())
        ]
        if malformed_rows:
            raise ValueError(f"Malformed TSV row in {path}")
        if require_nonempty and not path_rows:
            raise ValueError(f"Required replicate table is empty: {path}")
        expected_rep = None
        if expected_reps is not None:
            expected_rep = expected_reps[path_index]
            if "rep" not in schema or any(
                row["rep"] != str(expected_rep) for row in path_rows
            ):
                raise ValueError(f"replicate mismatch in {path}")
        if expected_chrom is not None:
            expected_chrom_value = str(expected_chrom)
            if "chrom" not in schema or any(
                row["chrom"] != expected_chrom_value for row in path_rows
            ):
                raise ValueError(f"Chromosome mismatch in {path}")
        validate_parquet_table(
            path.with_suffix(".parquet"),
            schema,
            len(path_rows),
            expected_rep=expected_rep,
            expected_chrom=expected_chrom
        )
        rows.extend(path_rows)
    if return_schema:
        return rows, expected_schema
    return rows


'''
internal: write a combined TSV and invoke the canonical Parquet writer.
'''
def _write_combined_stats_table(
    stats_dir,
    output_name,
    rows,
    fieldnames=None
):
    stats_path = Path(stats_dir)
    tsv_path = stats_path / f"{output_name}.tsv"
    parquet_path = stats_path / f"{output_name}.parquet"
    if fieldnames is None:
        fieldnames = list(rows[0]) if rows else []
    with open(tsv_path, "w", encoding="utf-8", newline="") as out_file:
        writer = csv.DictWriter(
            out_file,
            fieldnames=fieldnames,
            delimiter="\t",
            lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)
    converter = Path(__file__).parent / "write_parquet.py"
    subprocess.run(
        [sys.executable, str(converter), str(tsv_path), str(parquet_path)],
        check=True
    )


##### arguments ###############################################################
'''
define command-line arguments for
    - statistics output directory
    - number of replicates
    - list of chromosomes
    - chromosome index
    - list of tables to include / optional table selection
'''
parser = argparse.ArgumentParser()
parser.add_argument("--stats-dir", required=True)
parser.add_argument("--num-reps", type=int, required=True)
parser.add_argument("--chroms", nargs="+", required=True)
parser.add_argument("--chrom-index", type=int, required=True)
parser.add_argument(
    "--tables",
    nargs="+",
    choices=TABLE_CHOICES,
    default=TABLE_NAMES
)


##### main ####################################################################
if __name__ == "__main__":
    args = parser.parse_args()
    if args.chrom_index < 1 or args.chrom_index > len(args.chroms):
        raise ValueError(
            f"--chrom-index must be between 1 and {len(args.chroms)}"
        )

    stats_path = Path(args.stats_dir)
    chrom = args.chroms[args.chrom_index - 1]
    for table_name in args.tables:
        paths = [
            stats_path / f"{table_name}.rep_{rep}.chr{chrom}.tsv"
            for rep in range(1, args.num_reps + 1)
        ]
        try:
            rows, fieldnames = _combine_table_paths(
                paths,
                expected_reps=list(range(1, args.num_reps + 1)),
                expected_chrom=chrom,
                require_nonempty=table_name in NONEMPTY_TABLE_NAMES,
                return_schema=True
            )
        except FileNotFoundError as exc:
            raise FileNotFoundError(
                f"Missing chromosome files for {table_name}: {exc}"
            ) from exc
        _write_combined_stats_table(
            stats_path,
            f"{table_name}.chr{chrom}",
            rows,
            fieldnames=fieldnames
        )

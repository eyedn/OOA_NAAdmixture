###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           combine_sim_stats_tables.py
###############################################################################


from pathlib import Path
import csv
import pandas as pd


# read a TSV file into dictionary rows
def _read_tsv_rows(path):
    with open(path, "r", encoding="utf-8", newline="") as in_file:
        return list(csv.DictReader(in_file, delimiter="\t"))


# write dictionary rows to a TSV file
def _write_tsv_rows(path, rows):
    fieldnames = list(rows[0].keys()) if rows else []
    with open(path, "w", encoding="utf-8", newline="") as out_file:
        writer = csv.DictWriter(
            out_file,
            fieldnames=fieldnames,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


# write parquet output when pandas and pyarrow support are available
def _write_parquet(path, rows):
    if pd is None:
        return
    pd.DataFrame(rows).to_parquet(path, index=False)


# combine rows from required input paths and fail if any are missing
def _combine_paths(paths):
    missing_paths = [path for path in paths if not path.exists()]
    if missing_paths:
        missing = ", ".join(str(path) for path in missing_paths)
        raise FileNotFoundError(missing)

    rows = []
    for path in paths:
        rows.extend(_read_tsv_rows(path))
    return rows


# write one combined table with its requested output suffix
def _write_combined_table(stats_path, table_name, output_suffix, rows):
    out_prefix = f"{table_name}{output_suffix}"
    _write_tsv_rows(stats_path / f"{out_prefix}.tsv", rows)
    _write_parquet(stats_path / f"{out_prefix}.parquet", rows)


# combine per-replicate statistics tables after the stats array succeeds
def combine_sim_stats_tables(args, table_names):
    stats_path = Path(args.stats_dir)
    for table_name in table_names:
        per_rep_paths = [
            stats_path / f"{table_name}.rep_{rep}.tsv"
            for rep in range(1, args.num_reps + 1)
        ]
        try:
            rows = _combine_paths(per_rep_paths)
        except FileNotFoundError as exc:
            raise FileNotFoundError(
                f"Missing per-replicate files for {table_name}: {exc}"
            ) from exc
        _write_combined_table(stats_path, table_name, "", rows)


# combine per-replicate-per-chromosome statistics tables
def combine_chromosome_stats_tables(args, table_names):
    stats_path = Path(args.stats_dir)
    for table_name in table_names:
        chrom_paths = [
            stats_path / f"{table_name}.rep_{rep}.chr{chrom}.tsv"
            for rep in range(1, args.num_reps + 1)
            for chrom in args.chroms
        ]
        try:
            rows = _combine_paths(chrom_paths)
        except FileNotFoundError as exc:
            raise FileNotFoundError(
                f"Missing chromosome files for {table_name}: {exc}"
            ) from exc
        _write_combined_table(stats_path, table_name, ".chromosomes", rows)

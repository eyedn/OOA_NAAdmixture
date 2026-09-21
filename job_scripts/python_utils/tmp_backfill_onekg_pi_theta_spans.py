###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           tmp_backfill_onekg_pi_theta_spans.py
###############################################################################

# overview: rebuild and atomically republish empirical genome pi/theta spans.


##### set up ##################################################################
from pathlib import Path
import argparse
import csv
import math
import os
import shutil
import uuid
import pandas as pd
from onekg_utils.calc_onekg_stats import _aggregate_pi_theta_rows
from onekg_utils.read_tsv_rows import read_tsv_rows


TABLE_FAMILIES = (
    "pi_theta_stats_intergenic",
    "pi_theta_stats_full_callable_chrom"
)
GENOME_COLUMNS = [
    "rep", "pop", "stat", "value", "ne_value", "mutation_rate",
    "span", "segregating_sites", "wattersons_const"
]
NONSPAN_COLUMNS = [column for column in GENOME_COLUMNS if column != "span"]


##### internal functions #####################################################
'''
internal: require a nonempty table path before attempting to read its rows.
'''
def _require_table(path):
    if not path.is_file() or path.stat().st_size == 0:
        raise ValueError(f"Missing or empty table: {path}")


'''
internal: write one staged TSV/Parquet pair without modifying published paths.
'''
def _write_table_pair(prefix, rows):
    tsv_path = Path(f"{prefix}.tsv")
    parquet_path = Path(f"{prefix}.parquet")
    if not rows:
        raise ValueError(f"Cannot write an empty table: {prefix}")
    fieldnames = list(rows[0])
    if any(list(row) != fieldnames for row in rows):
        raise ValueError(f"Inconsistent row schema for {prefix}")
    with open(tsv_path, "w", encoding="utf-8", newline="") as out_file:
        writer = csv.DictWriter(
            out_file,
            fieldnames=fieldnames,
            delimiter="\t",
            lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)
    pd.DataFrame(rows, columns=fieldnames).to_parquet(parquet_path, index=False)


'''
internal: compare numeric values with strict tolerance and missing-value match.
'''
def _values_match(left, right):
    if pd.isna(left) and pd.isna(right):
        return True
    try:
        left_number = float(left)
        right_number = float(right)
    except (TypeError, ValueError):
        return left == right
    return math.isclose(
        left_number,
        right_number,
        rel_tol=1e-12,
        abs_tol=0.0
    )


'''
internal: check TSV and Parquet content agreement for an already-written pair.
'''
def _validate_tsv_parquet_agreement(prefix):
    tsv_path = Path(f"{prefix}.tsv")
    parquet_path = Path(f"{prefix}.parquet")
    _require_table(tsv_path)
    _require_table(parquet_path)
    tsv = pd.read_csv(tsv_path, sep="\t")
    parquet = pd.read_parquet(parquet_path)
    if tsv.columns.tolist() != parquet.columns.tolist():
        raise ValueError(f"TSV/Parquet schema mismatch for {prefix}")
    if len(tsv.index) != len(parquet.index):
        raise ValueError(f"TSV/Parquet row-count mismatch for {prefix}")
    for row_index in range(len(tsv.index)):
        for column in tsv.columns:
            if not _values_match(
                tsv.iloc[row_index][column],
                parquet.iloc[row_index][column]
            ):
                raise ValueError(
                    f"TSV/Parquet value mismatch for {prefix} "
                    f"row={row_index} column={column}"
                )
    return tsv


'''
internal: validate a genome table's schema, population/stat rows, and spans.
'''
def _validate_genome_table(prefix, expected_pops):
    table = _validate_tsv_parquet_agreement(prefix)
    if table.columns.tolist() != GENOME_COLUMNS:
        raise ValueError(f"Unexpected genome pi/theta schema for {prefix}")
    expected_keys = {
        (0, pop, stat)
        for pop in expected_pops
        for stat in ("pi", "theta")
    }
    actual_keys = {
        (int(row.rep), row.pop, row.stat)
        for row in table.itertuples(index=False)
    }
    if actual_keys != expected_keys or len(table.index) != len(expected_keys):
        raise ValueError(f"Unexpected genome pi/theta rows for {prefix}")
    spans = pd.to_numeric(table["span"], errors="coerce")
    if not spans.notna().all() or not spans.map(math.isfinite).all():
        raise ValueError(f"Genome pi/theta spans must be finite for {prefix}")
    if (spans <= 0).any():
        raise ValueError(f"Genome pi/theta spans must be positive for {prefix}")
    return table


'''
internal: validate one chromosome diversity source before genome aggregation.
'''
def _read_chromosome_rows(stats_dir, family, chrom, pop):
    path = stats_dir / f"{family}.rep_0.chr{chrom}.{pop}.tsv"
    _require_table(path)
    rows = read_tsv_rows(path)
    expected_keys = {
        ("0", str(chrom), pop, "pi"),
        ("0", str(chrom), pop, "theta")
    }
    actual_keys = {
        (str(row.get("rep")), str(row.get("chrom")), row.get("pop"),
         row.get("stat"))
        for row in rows
    }
    if actual_keys != expected_keys or len(rows) != len(expected_keys):
        raise ValueError(f"Unexpected chromosome pi/theta rows in {path}")
    return rows


'''
internal: rebuild one population's genome table from completed chromosome rows.
'''
def _build_population_rows(stats_dir, family, chroms, pop):
    chromosome_rows = []
    for chrom in chroms:
        chromosome_rows.extend(
            _read_chromosome_rows(stats_dir, family, chrom, pop)
        )
    return _aggregate_pi_theta_rows(chromosome_rows)


'''
internal: compare all non-span values with the completed published output.
'''
def _validate_nonspan_values(old_prefix, new_rows):
    old_table = _validate_tsv_parquet_agreement(old_prefix)
    unexpected_columns = set(old_table.columns) - set(GENOME_COLUMNS)
    if unexpected_columns:
        raise ValueError(f"Unexpected old genome schema for {old_prefix}")
    missing_columns = set(NONSPAN_COLUMNS) - set(old_table.columns)
    if missing_columns:
        raise ValueError(f"Incomplete old genome schema for {old_prefix}")
    old_by_key = {
        (int(row.rep), row.pop, row.stat): row
        for row in old_table.itertuples(index=False)
    }
    new_by_key = {
        (int(row["rep"]), row["pop"], row["stat"]): row
        for row in new_rows
    }
    if set(old_by_key) != set(new_by_key):
        raise ValueError(
            f"Old genome rows differ from rebuilt rows for {old_prefix}"
        )
    for key, old_row in old_by_key.items():
        new_row = new_by_key[key]
        for column in NONSPAN_COLUMNS:
            if not _values_match(getattr(old_row, column), new_row[column]):
                raise ValueError(
                    f"Non-span value changed for {old_prefix} key={key} "
                    f"column={column}"
                )


'''
internal: move staged files into place with backups and rollback on failure.
'''
def _promote_staged_pairs(stats_dir, staged_prefixes):
    backup_dir = stats_dir.parent / (
        f".{stats_dir.name}.pi_theta_span_backfill.backup.{uuid.uuid4().hex}"
    )
    backup_dir.mkdir()
    moved_backups = []
    promoted = []
    try:
        for prefix in staged_prefixes:
            for extension in ("tsv", "parquet"):
                target = stats_dir / f"{prefix.name}.{extension}"
                backup = backup_dir / target.name
                _require_table(target)
                os.replace(target, backup)
                moved_backups.append((target, backup))
        for prefix in staged_prefixes:
            for extension in ("tsv", "parquet"):
                staged = Path(f"{prefix}.{extension}")
                target = stats_dir / staged.name
                os.replace(staged, target)
                promoted.append(target)
    except Exception:
        for target in promoted:
            target.unlink(missing_ok=True)
        for target, backup in reversed(moved_backups):
            if backup.exists():
                os.replace(backup, target)
        _remove_temporary_dir(backup_dir)
        raise
    return backup_dir, moved_backups


'''
internal: discard temporary staging or backup directories created by this tool.
'''
def _remove_temporary_dir(path):
    if path.exists():
        shutil.rmtree(path)


##### main ####################################################################
'''
rebuild both empirical genome pi/theta families and republish validated pairs.
'''
def backfill_pi_theta_spans(stats_dir, chroms, pops):
    stats_dir = Path(stats_dir)
    if not stats_dir.is_dir():
        raise ValueError(f"Missing statistics directory: {stats_dir}")
    if not chroms or not pops:
        raise ValueError("At least one chromosome and population are required")
    stage_dir = stats_dir.parent / (
        f".{stats_dir.name}.pi_theta_span_backfill.stage.{uuid.uuid4().hex}"
    )
    stage_dir.mkdir()
    staged_prefixes = []
    try:
        for family in TABLE_FAMILIES:
            combined_rows = []
            for pop in pops:
                rows = _build_population_rows(stats_dir, family, chroms, pop)
                prefix_name = f"{family}.rep_0.{pop}"
                _validate_nonspan_values(stats_dir / prefix_name, rows)
                staged_prefix = stage_dir / prefix_name
                _write_table_pair(staged_prefix, rows)
                _validate_genome_table(staged_prefix, (pop,))
                staged_prefixes.append(staged_prefix)
                combined_rows.extend(rows)
            for prefix_name in (f"{family}.rep_0", family):
                _validate_nonspan_values(stats_dir / prefix_name, combined_rows)
                staged_prefix = stage_dir / prefix_name
                _write_table_pair(staged_prefix, combined_rows)
                _validate_genome_table(staged_prefix, pops)
                staged_prefixes.append(staged_prefix)

        backup_dir, moved_backups = _promote_staged_pairs(
            stats_dir,
            staged_prefixes
        )
        try:
            for staged_prefix in staged_prefixes:
                published_prefix = stats_dir / staged_prefix.name
                expected_pops = (
                    pops
                    if staged_prefix.name.endswith("rep_0")
                    or staged_prefix.name in TABLE_FAMILIES
                    else (staged_prefix.name.rsplit(".", 1)[1],)
                )
                _validate_genome_table(published_prefix, expected_pops)
        except Exception:
            for staged_prefix in staged_prefixes:
                for extension in ("tsv", "parquet"):
                    target = stats_dir / f"{staged_prefix.name}.{extension}"
                    target.unlink(missing_ok=True)
            for target, backup in reversed(moved_backups):
                if backup.exists():
                    os.replace(backup, target)
            _remove_temporary_dir(backup_dir)
            raise
        _remove_temporary_dir(backup_dir)
        return [stats_dir / prefix.name for prefix in staged_prefixes]
    finally:
        _remove_temporary_dir(stage_dir)


##### arguments ###############################################################
parser = argparse.ArgumentParser()
parser.add_argument("--stats-dir", required=True)
parser.add_argument("--chroms", nargs="+", required=True)
parser.add_argument("--pops", nargs="+", required=True)


if __name__ == "__main__":
    args = parser.parse_args()
    backfill_pi_theta_spans(args.stats_dir, args.chroms, args.pops)

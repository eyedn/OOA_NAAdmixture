###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           validate_parquet_table.py
###############################################################################

# overview: validate a replicate Parquet table against its canonical TSV.


##### set up ##################################################################
from pathlib import Path
import pandas as pd


##### main ####################################################################
def validate_parquet_table(
    parquet_path,
    tsv_schema,
    expected_row_count,
    expected_rep=None,
    expected_chrom=None,
):
    '''
    require a readable Parquet companion matching its TSV contract.
    '''
    parquet_path = Path(parquet_path)
    if not parquet_path.is_file() or parquet_path.stat().st_size == 0:
        raise ValueError(f"Missing or empty Parquet table: {parquet_path}")

    table = pd.read_parquet(parquet_path)
    if list(table.columns) != list(tsv_schema):
        raise ValueError(f"Parquet schema mismatch in {parquet_path}")
    if len(table.index) != expected_row_count:
        raise ValueError(f"Parquet row-count mismatch in {parquet_path}")

    for column, expected_value in (
        ("rep", expected_rep),
        ("chrom", expected_chrom),
    ):
        if expected_value is None:
            continue
        if column not in table.columns:
            raise ValueError(
                f"Parquet table lacks {column} in {parquet_path}"
            )
        observed_values = set(table[column].astype(str))
        if observed_values and observed_values != {str(expected_value)}:
            raise ValueError(
                f"Parquet {column} mismatch in {parquet_path}"
            )

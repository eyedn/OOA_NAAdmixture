###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           write_combined_stats_table.py
###############################################################################

# overview: write final combined statistics tables as TSV and Parquet.


##### set up ##################################################################
from pathlib import Path
import csv
import subprocess
import sys


##### main function ###########################################################
'''
write canonical TSV and Parquet outputs for one combined table. The TSV column
order comes from the first row and is reused by the Parquet conversion.
'''
def write_combined_stats_table(stats_dir, output_name, rows):
    stats_path = Path(stats_dir)
    tsv_path = stats_path / f"{output_name}.tsv"
    parquet_path = stats_path / f"{output_name}.parquet"
    fieldnames = list(rows[0]) if rows else []
    with open(tsv_path, "w", encoding="utf-8", newline="") as out_file:
        writer = csv.DictWriter(
            out_file,
            fieldnames=fieldnames,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)
    converter = Path(__file__).parents[1] / "write_parquet.py"
    subprocess.run(
        [sys.executable, str(converter), str(tsv_path), str(parquet_path)],
        check=True,
    )

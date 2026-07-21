###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           combine_table_paths.py
###############################################################################

# overview: read and concat. required TSVs without skipping missing inputs.


##### set up ##################################################################
from pathlib import Path
import csv


##### main function ###########################################################
'''
read and concatenate required TSV paths without silently skipping inputs.
Returns rows in the same order as the input path list.
'''
def combine_table_paths(paths):
    missing_paths = [Path(path) for path in paths if not Path(path).exists()]
    if missing_paths:
        raise FileNotFoundError(", ".join(str(path) for path in missing_paths))
    rows = []
    for path in paths:
        with open(path, "r", encoding="utf-8", newline="") as in_file:
            rows.extend(csv.DictReader(in_file, delimiter="\t"))
    return rows

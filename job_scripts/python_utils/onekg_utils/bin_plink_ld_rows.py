###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           bin_plink_ld_rows.py
###############################################################################
# bin PLINK pairwise LD rows into distance-decay summaries.


##### set up ##################################################################
from collections import defaultdict
import math


##### main function ###########################################################
'''
bin PLINK LD pairs into nonoverlapping genomic windows and distance bins.
Returns sufficient statistics and mean r2 for each retained bin.
'''
def bin_plink_ld_rows(
    pairs,
    rep,
    chrom,
    pop,
    window_size_bp=2_000_000,
    distance_bin_bp=5_000,
):
    grouped = defaultdict(lambda: {"sum_r2": 0.0, "n_pairs": 0})
    for pair in pairs:
        pos1_key = next(
            (key for key in ("POS_A", "POS1", "BP_A") if key in pair),
            None,
        )
        pos2_key = next(
            (key for key in ("POS_B", "POS2", "BP_B") if key in pair),
            None,
        )
        r2_key = next(
            (key for key in ("UNPHASED_R2", "R2") if key in pair),
            None,
        )
        if pos1_key is None or pos2_key is None or r2_key is None:
            raise ValueError("PLINK LD output lacks position or r2 columns")
        pos1 = int(pair[pos1_key])
        pos2 = int(pair[pos2_key])
        window_start = ((pos1 - 1) // window_size_bp) * window_size_bp
        if ((pos2 - 1) // window_size_bp) * window_size_bp != window_start:
            continue
        distance = pos2 - pos1
        if distance <= 0 or distance > window_size_bp:
            continue
        distance_bin = (
            ((distance - 1) // distance_bin_bp) + 1
        ) * distance_bin_bp
        key = (window_start, distance_bin)
        r2_value = float(pair[r2_key])
        grouped[key]["sum_r2"] += (
            r2_value if math.isfinite(r2_value) else 0.0
        )
        grouped[key]["n_pairs"] += 1
    output = []
    for (window_start, distance_bin), values in sorted(grouped.items()):
        output.append(
            {
                "rep": rep,
                "chrom": chrom,
                "pop": pop,
                "window_start": window_start,
                "window_end": window_start + window_size_bp,
                "distance_bin_bp": distance_bin,
                "mean_r2": values["sum_r2"] / values["n_pairs"],
                "sum_r2": values["sum_r2"],
                "n_pairs": values["n_pairs"],
            }
        )
    return output

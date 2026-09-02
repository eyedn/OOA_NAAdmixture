###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           select_simOnekgDownsample_snps.py
###############################################################################

# overview: select simulation SNPs and write their empirical-density contract.


##### set up ##################################################################
from pathlib import Path
import argparse
import json
import sys

from shared_utils.write_stats_table import write_stats_table
from simOnekgDownsample_utils import (
    parse_sim_positions,
    parse_vcftools_density,
    select_snp_density
)


##### arguments ###############################################################
parser = argparse.ArgumentParser()
parser.add_argument("--empirical-density", required=True)
parser.add_argument("--simulation-density", required=True)
parser.add_argument("--simulation-positions", required=True)
parser.add_argument("--wanted-snps", required=True)
parser.add_argument("--contract-json", required=True)
parser.add_argument("--stats-dir", required=True)
parser.add_argument("--window-size-bp", required=True, type=int)
parser.add_argument("--seed", required=True, type=int)
parser.add_argument("--rep", required=True, type=int)
parser.add_argument("--chrom", required=True)


##### main ####################################################################
if __name__ == "__main__":
    args = parser.parse_args()
    with open(args.empirical_density, "r", encoding="utf-8") as in_file:
        empirical_rows = parse_vcftools_density(in_file)
    with open(args.simulation_density, "r", encoding="utf-8") as in_file:
        simulation_rows = parse_vcftools_density(in_file)
    with open(args.simulation_positions, "r", encoding="utf-8") as in_file:
        simulation_positions = parse_sim_positions(in_file)
    selected, contracted_rows = select_snp_density(
        empirical_rows,
        simulation_rows,
        simulation_positions,
        args.window_size_bp,
        args.seed,
        args.rep
    )
    contracted_chroms = {str(row["chrom"]) for row in contracted_rows}
    expected_chrom = str(args.chrom)
    if expected_chrom.lower().startswith("chr"):
        expected_chrom = expected_chrom[3:]
    if contracted_chroms != {expected_chrom}:
        raise ValueError("Density contract chromosome mismatch")

    with open(args.wanted_snps, "w", encoding="utf-8") as out_file:
        out_file.writelines(
            f"{chrom}\t{position}\n" for chrom, position in selected
        )
    with open(args.contract_json, "w", encoding="utf-8") as out_file:
        json.dump(contracted_rows, out_file, indent=2)
        out_file.write("\n")
    stats_dir = Path(args.stats_dir)
    stats_dir.mkdir(parents=True, exist_ok=True)
    write_stats_table(
        stats_dir / f"snp_density.rep_{args.rep}.chr{args.chrom}",
        contracted_rows
    )
    for row in contracted_rows:
        if row["simulation_below_target"]:
            print(
                "WARNING: simulation SNP deficit "
                f"rep={row['rep']} chr={row['chrom']} "
                f"bin_start={row['bin_start']} "
                f"deficit={row['deficit_snp_count']}",
                file=sys.stderr
            )

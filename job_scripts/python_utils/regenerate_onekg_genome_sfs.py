###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           regenerate_onekg_genome_sfs.py
###############################################################################

import argparse
from pathlib import Path

from onekg_utils.calc_onekg_stats import _aggregate_projected_sfs_rows
from onekg_utils.read_tsv_rows import read_tsv_rows
from shared_utils.write_stats_table import write_stats_table


parser = argparse.ArgumentParser()
parser.add_argument("--stats-dir", required=True)
parser.add_argument("--chroms", nargs="+", required=True)
parser.add_argument("--pops", nargs="+", required=True)


if __name__ == "__main__":
    args = parser.parse_args()
    stats_dir = Path(args.stats_dir)
    genome_by_pop = {}
    for pop in args.pops:
        chromosome_rows = []
        for chrom in args.chroms:
            chromosome_rows.extend(
                read_tsv_rows(
                    stats_dir / f"sfs.rep_0.chr{chrom}.{pop}.tsv"
                )
            )
        genome_by_pop[pop] = _aggregate_projected_sfs_rows(chromosome_rows)
    for pop, rows in genome_by_pop.items():
        write_stats_table(stats_dir / f"sfs.rep_0.{pop}", rows)
    combined = [
        row
        for pop in args.pops
        for row in genome_by_pop[pop]
    ]
    write_stats_table(stats_dir / "sfs.rep_0", combined)
    write_stats_table(stats_dir / "sfs", combined)

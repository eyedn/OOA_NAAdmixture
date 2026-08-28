###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           regenerate_onekg_chr_sfs.py
###############################################################################

import argparse
import gzip
from pathlib import Path

from onekg_utils.calc_onekg_stats import _build_projected_folded_1d_sfs_rows
from onekg_utils.read_onekg_sample_pops import read_onekg_sample_pops
from onekg_utils.scan_onekg_vcf import scan_onekg_vcf
from shared_utils.write_stats_table import write_stats_table


parser = argparse.ArgumentParser()
parser.add_argument("--vcf-path", required=True)
parser.add_argument("--unrels-path", required=True)
parser.add_argument("--fam-path", required=True)
parser.add_argument("--stats-dir", required=True)
parser.add_argument("--chrom", required=True)
parser.add_argument("--sfs-size", type=int, required=True)
parser.add_argument("--sfs-size-pop-ref", required=True)
parser.add_argument("--pops", nargs="+", required=True)


if __name__ == "__main__":
    args = parser.parse_args()
    sample_pops = read_onekg_sample_pops(
        args.unrels_path, args.fam_path, args.pops
    )
    open_vcf = gzip.open if args.vcf_path.endswith(".gz") else open
    with open_vcf(args.vcf_path, "rt", encoding="utf-8") as vcf_file:
        scan = scan_onekg_vcf(vcf_file, sample_pops, args.pops)
    stats_dir = Path(args.stats_dir)
    combined = []
    by_pop = {}
    for pop in args.pops:
        rows = _build_projected_folded_1d_sfs_rows(
            0,
            args.chrom,
            pop,
            list(scan["counts_by_pop"][pop].values()),
            scan["sample_counts"][pop],
            2 * args.sfs_size,
            args.sfs_size_pop_ref
        )
        by_pop[pop] = rows
        combined.extend(rows)
    for pop, rows in by_pop.items():
        write_stats_table(
            stats_dir / f"sfs.rep_0.chr{args.chrom}.{pop}", rows
        )
    write_stats_table(stats_dir / f"sfs.rep_0.chr{args.chrom}", combined)
    write_stats_table(stats_dir / f"sfs.chr{args.chrom}", combined)

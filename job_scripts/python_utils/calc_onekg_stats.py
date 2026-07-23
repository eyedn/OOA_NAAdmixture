###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           calc_onekg_stats.py
###############################################################################

# overview: calculate 1000 Genomes chromosome or genome statistics from CLI 
# inputs.


##### set up ##################################################################
import argparse
from onekg_utils.calc_onekg_stats import calc_onekg_stats


##### arguments ###############################################################
'''
define shared KING and output arguments plus the inputs required only for
chromosome-level or genome-level empirical statistics.
'''
parser = argparse.ArgumentParser()
parser.add_argument(
    "--analysis-level",
    choices=["chromosome", "genome"],
    required=True,
)
parser.add_argument("--stats-dir", required=True)
parser.add_argument("--pop", required=True)
parser.add_argument("--king-path", required=True)
parser.add_argument("--vcf-path")
parser.add_argument("--ld-vcf-path")
parser.add_argument("--intergenic-vcf-path")
parser.add_argument("--intergenic-bed-path")
parser.add_argument("--span-incl-bed-path")
parser.add_argument("--unrels-path")
parser.add_argument("--fam-path")
parser.add_argument("--chr-lens-path")
parser.add_argument("--chrom")
parser.add_argument("--chroms", nargs="+")
parser.add_argument("--mutation-rate", type=float)
parser.add_argument("--ld-decay-window-size-bp", type=int)
parser.add_argument("--ld-decay-distance-bin-bp", type=int)
parser.add_argument("--ld-decay-maf-threshold", type=float)
parser.add_argument("--pops", nargs="+")


##### main ####################################################################
if __name__ == "__main__":
    args = parser.parse_args()
    calc_onekg_stats(args)

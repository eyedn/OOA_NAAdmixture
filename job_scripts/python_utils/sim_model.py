###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           sim_model.py
###############################################################################

# overview: parse simulation request and hand it to the model runner.


##### set up ##################################################################
import argparse
from sim_utils.run_simulation import run_simulation


##### arguments ###############################################################
'''
define command-line arguments for
    - tree, vcf, metadata, ancestry paths
    - samples size and random seed
    - chromosome and genetic map
    - msprime simulation model
    - model parameters and labels
'''
parser = argparse.ArgumentParser()
parser.add_argument("--tree-prefix", required=True)
parser.add_argument("--pickle-prefix", required=True)
parser.add_argument("--vcf-path", required=True)
parser.add_argument("--sample-metadata-path", required=True)
parser.add_argument("--anc-dir", required=True)
parser.add_argument("--global-anc-dir", required=True)
parser.add_argument("--sample-size", type=int, required=True)
parser.add_argument("--seed", type=int, required=True)
parser.add_argument("--msprime-model", required=True)
parser.add_argument("--chromosome", required=True)
parser.add_argument("--genetic-map", required=True)
parser.add_argument("--generation-time", type=float, required=True)
parser.add_argument("--mutation-rate", type=float, required=True)
parser.add_argument("--t-af-years", type=float, required=True)
parser.add_argument("--t-ooa-years", type=float, required=True)
parser.add_argument("--t-eu0-years", type=float, required=True)
parser.add_argument("--t-eg-years", type=float, required=True)
parser.add_argument("--r-eu0", type=float, required=True)
parser.add_argument("--r-eu", type=float, required=True)
parser.add_argument("--r-af", type=float, required=True)
parser.add_argument("--n-a", type=float, required=True)
parser.add_argument("--n-af1", type=float, required=True)
parser.add_argument("--n-b", type=float, required=True)
parser.add_argument("--n-eu0", type=float, required=True)
parser.add_argument("--m-af-b", type=float, required=True)
parser.add_argument("--m-af-eu", type=float, required=True)
parser.add_argument("--admixture-time", type=float, required=True)
parser.add_argument("--admix-generation-count", type=int, required=True)
parser.add_argument("--admix-mixing-generation-count", type=int, required=True)
parser.add_argument("--admix-ne-by-generation", required=True)
parser.add_argument("--admix-afr-props-by-generation", required=True)
parser.add_argument("--admix-eur-props-by-generation", required=True)
parser.add_argument("--admix-prioradmix-props-by-generation", required=True)
parser.add_argument("--admix-modern-growth-rate", type=float, required=True)
parser.add_argument("--census-time-offset", type=float, required=True)
parser.add_argument("--pops", nargs="+", required=True)


##### main ####################################################################
if __name__ == "__main__":
    args = parser.parse_args()
    run_simulation(args)

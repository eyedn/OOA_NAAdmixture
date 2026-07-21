###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           write_pop_subset.py
###############################################################################
# write one population's deterministic PLINK and KING keep-file subset.


##### set up ##################################################################
import argparse


##### arguments ###############################################################
'''
define command-line arguments for
    - pop to subset
    - sample size and list of all pops
'''
parser = argparse.ArgumentParser()
parser.add_argument("--subset-path", required=True)
parser.add_argument("--pop", required=True)
parser.add_argument("--sample-size", type=int, required=True)
parser.add_argument("--pops", nargs="+", required=True)


##### main ####################################################################
if __name__ == "__main__":
    args = parser.parse_args()
    subset_path = args.subset_path
    pop = args.pop
    pops = args.pops
    sample_size = args.sample_size

    # find the target population's sample indexes in the full population order.
    pop_idx = pops.index(pop)
    start = pop_idx * sample_size
    end = (pop_idx + 1) * sample_size

    '''
    write one PLINK-style KING keep-file row per sample, using:
        - family ID 0
        - sample IDs formatted as <population>_1, <population>_2, and so on
    '''
    with open(subset_path, "w", encoding="utf-8") as out_file:
        for sample_idx in range(start, end):
            sample_name = f"{pop}_{sample_idx - start + 1}"
            out_file.write(f"0\t{sample_name}\n")

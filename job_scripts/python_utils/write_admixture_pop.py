###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           write_admixture_pop.py
###############################################################################
# write supervised labels in the exact final FAM sample order.


##### set up ##################################################################
import argparse
from stats_utils.read_fam_order import read_fam_order
from stats_utils.read_sample_metadata import read_sample_metadata


##### arguments ###############################################################
'''
define command-line arguments for
    - sample/population metadata from sim
    - path to read plink2 .fam file
    - path to output .pop file for supervised ADMIXTURE
'''
parser = argparse.ArgumentParser()
parser.add_argument("--sample-metadata-path", required=True)
parser.add_argument("--fam-path", required=True)
parser.add_argument("--pop-path", required=True)


##### main ####################################################################
if __name__ == "__main__":
    args = parser.parse_args()
    sample_metadata_path = args.sample_metadata_path
    fam_path = args.fam_path
    pop_path = args.pop_path

    # read sample/population metadata from simulation with read_sample_metadata
    metadata_by_sample = {
        (row["fid"], row["iid"]): row
        for row in read_sample_metadata(sample_metadata_path)
    }

    '''
    - reads the final .fam file in its actual sample order
    - writes one pop label per line to the .pop file,
    preserving the exact sample order expected by ADMIXTURE
    '''
    with open(pop_path, "w", encoding="utf-8") as out_file:
        for row in read_fam_order(fam_path):
            sample_key = (row["fid"], row["iid"])
            if sample_key not in metadata_by_sample:
                raise ValueError(
                    f"Missing sample metadata for: {row['iid']}"
                )
            label = metadata_by_sample[sample_key]["supervised_label"]
            out_file.write(f"{label}\n")

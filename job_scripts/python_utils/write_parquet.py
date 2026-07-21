###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           write_parquet.py
###############################################################################

# overview: convert canonical TSV statistics tables to Parquet.


##### set up ##################################################################
import argparse
import pandas as pd


##### arguments ###############################################################
'''
define the source TSV and destination Parquet paths.
'''
parser = argparse.ArgumentParser()
parser.add_argument("input_tsv")
parser.add_argument("output_parquet")


##### main ####################################################################
if __name__ == "__main__":
    args = parser.parse_args()
    table = pd.read_csv(args.input_tsv, sep="\t")
    table.to_parquet(args.output_parquet, index=False)

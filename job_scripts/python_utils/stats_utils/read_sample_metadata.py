###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           read_sample_metadata.py
###############################################################################


# pattern: Imperative Shell


def read_sample_metadata(sample_metadata_path):
    import csv

    with open(sample_metadata_path, "r", encoding="utf-8") as in_file:
        return list(csv.DictReader(in_file, delimiter="\t"))

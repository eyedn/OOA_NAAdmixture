###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           subset_pairs.py
###############################################################################


def subset_pairs(subset_path, pop, pops, sample_size):
    pop_idx = pops.index(pop)
    start = pop_idx * sample_size
    end = (pop_idx + 1) * sample_size
    with open(subset_path, "w", encoding="utf-8") as out_file:
        for sample_idx in range(start, end):
            sample_name = f"{pop}_{sample_idx - start + 1}"
            out_file.write(f"0\t{sample_name}\n")

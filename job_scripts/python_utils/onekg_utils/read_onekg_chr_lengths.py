###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           read_onekg_chr_lengths.py
###############################################################################
# read 1000 Genomes chromosome lengths for statistics normalization.


##### main function ###########################################################
'''
read full and post-QC sequence lengths for one chromosome. Returns the callable
length followed by the full chromosome length used by diversity summaries.
'''
def read_onekg_chr_lengths(path, chrom):
    with open(path, "r", encoding="utf-8") as in_file:
        rows = [line.split() for line in in_file if line.strip()]
    if not rows:
        raise ValueError("Chromosome-length file is empty")
    headers = [value.lower() for value in rows[0]]
    chrom_names = ("chrom", "chr", "chromosome")
    full_names = ("chr_len", "length", "full_chrom_length")
    qc_names = ("chr_len_after_qc", "length_after_qc", "callable_length")
    chrom_index = next(
        (headers.index(name) for name in chrom_names if name in headers),
        None,
    )
    full_index = next(
        (headers.index(name) for name in full_names if name in headers),
        None,
    )
    qc_index = next(
        (headers.index(name) for name in qc_names if name in headers),
        None,
    )
    if None in (chrom_index, full_index, qc_index):
        raise ValueError(
            "Chromosome-length file needs chrom, chr_len, and "
            "chr_len_after_qc columns"
        )
    target = str(chrom).removeprefix("chr")
    for fields in rows[1:]:
        row_chrom = fields[chrom_index].removeprefix("chr")
        if row_chrom == target:
            return float(fields[qc_index]), float(fields[full_index])
    raise ValueError(f"Chromosome {chrom} is absent from chromosome lengths")

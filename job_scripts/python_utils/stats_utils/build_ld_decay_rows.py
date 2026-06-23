###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           build_ld_decay_rows.py
###############################################################################


import numpy as np
import allel


# calculate ld decay and return results per population
def build_ld_decay_rows(
    ts,
    rep,
    pops,
    sample_size,
    window_size_bp=2_000_000,
    distance_bin_bp=5_000,
    maf_threshold=0.10,
):
    rows = []
    positions = np.array([site.position for site in ts.sites()])

    # get genotype matrix from tree sequence
    haploid_genotypes = ts.genotype_matrix().astype(int)
    if haploid_genotypes.shape[1] % 2 != 0:
        raise ValueError("Expected diploid genotypes with an even node count")

    # reshape resulting matrix for diploid genotypes
    num_individuals = haploid_genotypes.shape[1] // 2
    genotypes = haploid_genotypes.reshape(
        haploid_genotypes.shape[0],
        num_individuals,
        2,
    )
    requested_individuals = len(pops) * sample_size
    if requested_individuals > num_individuals:
        raise ValueError(
            "Requested populations require more individuals than available"
        )

    # fitler snps by MAF
    genotype_array = allel.GenotypeArray(genotypes[:, :requested_individuals])
    allele_counts = genotype_array.count_alleles()
    alt_freq = allele_counts[:, 1] / allele_counts.sum(axis=1)
    maf = np.minimum(alt_freq, 1 - alt_freq)
    keep = maf >= maf_threshold
    positions = positions[keep]

    # convert matrix to 0/1/2 counts for sites that satisfied the maf threshold
    alt_counts = genotype_array[keep].to_n_alt()

    for pop in pops:
        # subset pop-specific individuals
        pop_idx = pops.index(pop)
        start = pop_idx * sample_size
        end = start + sample_size
        pop_alt_counts = alt_counts[:, start:end]

        # divide the contig into non-overlapping 2Mb windows
        for window_start in range(0, int(ts.sequence_length), window_size_bp):
            window_end = min(window_start + window_size_bp, ts.sequence_length)
            in_window = np.where(
                (positions >= window_start) & (positions < window_end)
            )[0]
            if len(in_window) < 2:
                continue

            # calculate r and r2 for all pairs in this 2Mb window
            window_positions = positions[in_window]
            block = pop_alt_counts[in_window, :]
            r_values = allel.rogers_huff_r(block)
            row_idx, col_idx = np.triu_indices(len(window_positions), k=1)
            distances = window_positions[col_idx] - window_positions[row_idx]
            r2_values = np.asarray(r_values) ** 2
            valid = (distances > 0) & (distances <= window_size_bp)
            distances = distances[valid]
            r2_values = r2_values[valid]

            # bin all mutation pairs by their distance and summarize r2
            bin_values = {}
            distance_bins = (
                ((distances.astype(int) - 1) // distance_bin_bp) + 1
            ) * distance_bin_bp
            for distance_bin, r2_value in zip(distance_bins, r2_values):
                bin_values.setdefault(int(distance_bin), []).append(
                    float(r2_value)
                )

            for distance_bin, r2_values in sorted(bin_values.items()):
                rows.append(
                    {
                        "rep": rep,
                        "pop": pop,
                        "window_start": int(window_start),
                        "window_end": int(window_end),
                        "distance_bin_bp": int(distance_bin),
                        "mean_r2": float(np.nanmean(r2_values)),
                        "sum_r2": float(np.nansum(r2_values)),
                        "n_pairs": len(r2_values),
                    }
                )
    return rows

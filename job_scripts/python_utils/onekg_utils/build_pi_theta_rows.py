###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           build_pi_theta_rows.py
###############################################################################

# overview: build pop. pi and theta rows from scanned VCF genotype counts.


##### main function ###########################################################
'''
build pi and Watterson theta rows using one explicit sequence denominator.
Returns separate rows because the statistics use different estimators.
'''
def build_pi_theta_rows(
    rep,
    chrom,
    pop,
    allele_counts,
    span,
    mutation_rate,
):
    if span <= 0:
        raise ValueError("Sequence span must be positive")
    sample_sizes = {ref_count + alt_count
                    for ref_count, alt_count in allele_counts}
    if len(sample_sizes) != 1:
        raise ValueError("Allele counts must have one called sample size")
    num_haplotypes = sample_sizes.pop()
    if num_haplotypes < 2:
        raise ValueError("At least two called haplotypes are required")

    pi_numerator = sum(
        2 * ref_count * alt_count
        / (num_haplotypes * (num_haplotypes - 1))
        for ref_count, alt_count in allele_counts
    )
    segregating_site_count = sum(
        ref_count > 0 and alt_count > 0
        for ref_count, alt_count in allele_counts
    )
    wattersons_const = sum(
        1 / value for value in range(1, num_haplotypes)
    )
    pi_value = pi_numerator / span
    segregating_sites = segregating_site_count / span
    theta_value = segregating_sites / wattersons_const
    common = {
        "rep": rep,
        "chrom": chrom,
        "pop": pop,
        "mutation_rate": mutation_rate,
        "span": span,
        "wattersons_const": wattersons_const,
    }
    return [
        {
            **common,
            "stat": "pi",
            "value": pi_value,
            "ne_value": pi_value / (4 * mutation_rate),
            "segregating_sites": None,
        },
        {
            **common,
            "stat": "theta",
            "value": theta_value,
            "ne_value": theta_value / (4 * mutation_rate),
            "segregating_sites": segregating_sites,
        },
    ]

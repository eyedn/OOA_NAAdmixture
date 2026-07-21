###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           aggregate_pi_theta_rows.py
###############################################################################

# overview: aggregate chromosome pi and theta summaries at genome scope.


##### set up ##################################################################
from collections import defaultdict


##### main function ###########################################################
'''
aggregate chromosome pi and theta rows at genome scope. Pi receives its
sample-derived finite-sample correction, while theta uses segregating sites.
'''
def aggregate_pi_theta_rows(rows, haplotypes_by_pop):
    grouped = defaultdict(list)
    for row in rows:
        grouped[(int(row["rep"]), row["pop"], row["stat"])].append(row)

    aggregated = []
    for (rep, pop, stat), stat_rows in sorted(grouped.items()):
        span = sum(float(row["span"]) for row in stat_rows)
        mutation_rate = float(stat_rows[0]["mutation_rate"])
        wattersons_const = float(stat_rows[0]["wattersons_const"])
        if stat == "pi":
            value = sum(
                float(row["value"]) * float(row["span"])
                for row in stat_rows
            ) / span
            num_haplotypes = haplotypes_by_pop.get(pop, 0)
            if num_haplotypes < 2:
                raise ValueError(f"Invalid haplotype count for pop={pop}")
            value *= num_haplotypes / (num_haplotypes - 1)
            segregating_sites = None
        elif stat == "theta":
            segregating_sites = sum(
                float(row["segregating_sites"]) * float(row["span"])
                for row in stat_rows
            ) / span
            value = segregating_sites / wattersons_const
        else:
            raise ValueError(f"Unknown statistic {stat}")
        aggregated.append(
            {
                "rep": rep,
                "pop": pop,
                "stat": stat,
                "value": value,
                "ne_value": value / (4 * mutation_rate),
                "mutation_rate": mutation_rate,
                "segregating_sites": segregating_sites,
                "wattersons_const": wattersons_const,
            }
        )
    return aggregated

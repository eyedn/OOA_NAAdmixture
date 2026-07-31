###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           aggregate_genome_stats.py
###############################################################################

# overview: aggregate chromosome stats. and related handoffs at genome scope.


##### set up ##################################################################
from collections import defaultdict
from pathlib import Path
import csv
import math
import pandas as pd
from shared_utils import log_msg
from .parse_king_file import parse_king_file
from .read_fam_order import read_q_rows


##### internal functions ######################################################
'''
return whether a table value can be included in numeric summaries.
'''
def _finite(value):
    return value is not None and not math.isnan(float(value))


'''
read one TSV as dictionaries while preserving column names and row order.
'''
def _read_tsv_rows(path):
    with open(path, "r", encoding="utf-8", newline="") as in_file:
        return list(csv.DictReader(in_file, delimiter="\t"))


'''
write dictionary rows to a TSV using the first row's column order.
'''
def _write_tsv_rows(path, rows):
    fieldnames = list(rows[0].keys()) if rows else []
    with open(path, "w", encoding="utf-8", newline="") as out_file:
        writer = csv.DictWriter(
            out_file,
            fieldnames=fieldnames,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


'''
convert pandas or simple table objects to a list of row dictionaries.
'''
def _table_to_rows(table):
    if hasattr(table, "to_dict"):
        return table.to_dict("records")
    return list(table)


'''
read one table type across required chromosome-level outputs in input order.
'''
def _read_chrom_tables(stats_dir, table_name, rep, chroms):
    rows = []
    for chrom in chroms:
        path = stats_dir / f"{table_name}.rep_{rep}.chr{chrom}.tsv"
        log_msg(f"reading chromosome table path={path}")
        rows.extend(_read_tsv_rows(path))
    return rows


'''
write one genome-level table as TSV and, when available, Parquet.
'''
def _write_table(stats_dir, table_name, rep, rows):
    out_tsv = stats_dir / f"{table_name}.rep_{rep}.tsv"
    out_parquet = stats_dir / f"{table_name}.rep_{rep}.parquet"
    log_msg(f"writing genome table path={out_tsv}")
    _write_tsv_rows(out_tsv, rows)
    if pd is not None:
        pd.DataFrame(rows).to_parquet(out_parquet, index=False)


'''
aggregate chromosome ancestry estimates using represented sequence span as the
weight. Returns one genome-level row per sample.
'''
def _aggregate_ancestry_rows(rows):
    grouped = {}
    for row in rows:
        key = (
            int(row["rep"]),
            row["pop"],
            row["sample_id"],
            row["vcf_sample_id"],
        )
        current = grouped.setdefault(
            key,
            {
                "rep": int(row["rep"]),
                "pop": row["pop"],
                "sample_id": row["sample_id"],
                "vcf_sample_id": row["vcf_sample_id"],
                "afr_tspop_weighted": 0.0,
                "eur_tspop_weighted": 0.0,
                "span": 0.0,
            },
        )
        span = float(row.get("span", 1.0))
        current["afr_tspop_weighted"] += float(row["afr_tspop"]) * span
        current["eur_tspop_weighted"] += float(row["eur_tspop"]) * span
        current["span"] += span

    aggregated = []
    for row in grouped.values():
        span = row["span"]
        aggregated.append(
            {
                "rep": row["rep"],
                "pop": row["pop"],
                "sample_id": row["sample_id"],
                "vcf_sample_id": row["vcf_sample_id"],
                "afr_tspop": row["afr_tspop_weighted"] / span,
                "eur_tspop": row["eur_tspop_weighted"] / span,
                "span": span,
            }
        )
    return sorted(
        aggregated,
        key=lambda row: (row["rep"], row["pop"], row["sample_id"]),
    )


'''
sum one-dimensional SFS bins across chromosomes for each population.
'''
def _aggregate_sfs_rows(rows):
    grouped = defaultdict(float)
    for row in rows:
        key = (
            int(row["rep"]),
            row["pop"],
            int(row["derived_allele_count"]),
        )
        grouped[key] += float(row["count"])

    return [
        {
            "rep": rep,
            "pop": pop,
            "derived_allele_count": derived_allele_count,
            "count": count,
        }
        for (rep, pop, derived_allele_count), count in sorted(grouped.items())
    ]


'''
sum pairwise two-dimensional SFS bins across chromosomes.
'''
def _aggregate_sfs_2d_rows(rows):
    grouped = defaultdict(float)
    for row in rows:
        key = (
            int(row["rep"]),
            row["pop1"],
            row["pop2"],
            int(row["pop1_count"]),
            int(row["pop2_count"]),
        )
        grouped[key] += float(row["count"])

    return [
        {
            "rep": rep,
            "pop1": pop1,
            "pop2": pop2,
            "pop1_count": pop1_count,
            "pop2_count": pop2_count,
            "count": count,
        }
        for (
            rep,
            pop1,
            pop2,
            pop1_count,
            pop2_count,
        ), count in sorted(grouped.items())
    ]


'''
aggregate pi by sequence span and theta by segregating sites. Pi receives the
configured finite-sample correction from all sampled haploid chromosomes.
'''
def _aggregate_pi_theta_rows(rows, sample_size):
    num_haplotypes = 2 * sample_size
    if num_haplotypes < 2:
        raise ValueError("Sample size must include at least one individual")

    grouped = defaultdict(list)
    for row in rows:
        grouped[(int(row["rep"]), row["pop"], row["stat"])].append(row)

    aggregated = []
    for (rep, pop, stat), stat_rows in sorted(grouped.items()):
        mutation_rate = float(stat_rows[0]["mutation_rate"])
        total_span = sum(float(row["span"]) for row in stat_rows)
        if stat == "pi":
            value = sum(
                float(row["value"]) * float(row["span"])
                for row in stat_rows
            ) / total_span
            value *= num_haplotypes / (num_haplotypes - 1)
            segregating_sites = None
            wattersons_const = None
        elif stat == "theta":
            segregating_sites = sum(
                float(row["segregating_sites"]) * float(row["span"])
                for row in stat_rows
                if _finite(row["segregating_sites"])
            ) / total_span
            wattersons_const = float(stat_rows[0]["wattersons_const"])
            value = segregating_sites / wattersons_const
        else:
            raise ValueError(f"Unknown pi/theta stat: {stat}")

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


'''
combine LD-decay bins across chromosome windows using summed r2 and pair counts.
'''
def _aggregate_ld_decay_rows(rows):
    grouped = defaultdict(lambda: {"sum_r2": 0.0, "n_pairs": 0})
    for row in rows:
        key = (int(row["rep"]), row["pop"], int(row["distance_bin_bp"]))
        grouped[key]["sum_r2"] += float(row["sum_r2"])
        grouped[key]["n_pairs"] += int(row["n_pairs"])

    aggregated = []
    for (rep, pop, distance_bin_bp), values in sorted(grouped.items()):
        n_pairs = values["n_pairs"]
        mean_r2 = values["sum_r2"] / n_pairs if n_pairs else math.nan
        aggregated.append(
            {
                "rep": rep,
                "pop": pop,
                "distance_bin_bp": distance_bin_bp,
                "mean_r2": mean_r2,
                "sum_r2": values["sum_r2"],
                "n_pairs": n_pairs,
            }
        )
    return aggregated


##### main function ###########################################################
'''
aggregate chromosome statistics and genome-level KING and ADMIXTURE handoffs
into one table set for a simulation replicate.
'''
def aggregate_genome_stats(args):
    stats_dir = Path(args.stats_dir)
    admixture_dir = Path(args.admixture_dir)
    king_dir = Path(args.king_dir)

    # aggregate span-weighted ancestry, then join ADMIXTURE in final FAM order.
    log_msg(f"aggregating chromosome stats rep={args.rep}")
    ancestry_rows = _read_chrom_tables(
        stats_dir,
        "ancestry",
        args.rep,
        args.chroms,
    )
    ancestry = _aggregate_ancestry_rows(ancestry_rows)

    genome_prefix = f"{args.genetic_map}_{args.rep}_genome_all"
    q_path = admixture_dir / f"{genome_prefix}.supervised.2.Q"
    fam_path = admixture_dir / f"{genome_prefix}.fam"
    q_rows = read_q_rows(args.rep, q_path, fam_path)
    q_by_sample = {
        (row["pop"], row["vcf_sample_id"]): row
        for row in q_rows
    }
    for row in ancestry:
        q_row = q_by_sample.get((row["pop"], row["vcf_sample_id"]), {})
        row["afr_q"] = q_row.get("afr_q")
        row["eur_q"] = q_row.get("eur_q")
    _write_table(stats_dir, "ancestry", args.rep, ancestry)

    # combine population-specific genome KING outputs into one kinship table.
    king_rows = []
    for pop in args.pops:
        king_path = (
            king_dir / f"{args.genetic_map}_{args.rep}_genome_{pop}.kin0"
        )
        king_table = parse_king_file(king_path, args.rep, pop)
        king_rows.extend(_table_to_rows(king_table))
    _write_table(stats_dir, "kinship", args.rep, king_rows)

    # aggregate chromosome summaries using each statistic's preserved contract.
    pi_theta = _read_chrom_tables(
        stats_dir,
        "pi_theta_stats",
        args.rep,
        args.chroms,
    )
    _write_table(
        stats_dir,
        "pi_theta_stats",
        args.rep,
        _aggregate_pi_theta_rows(pi_theta, args.sample_size),
    )

    sfs = _read_chrom_tables(stats_dir, "sfs", args.rep, args.chroms)
    _write_table(stats_dir, "sfs", args.rep, _aggregate_sfs_rows(sfs))

    sfs_2d = _read_chrom_tables(stats_dir, "sfs_2d", args.rep, args.chroms)
    _write_table(
        stats_dir,
        "sfs_2d",
        args.rep,
        _aggregate_sfs_2d_rows(sfs_2d),
    )

    ld_decay = _read_chrom_tables(stats_dir, "ld_decay", args.rep, args.chroms)
    _write_table(
        stats_dir,
        "ld_decay",
        args.rep,
        _aggregate_ld_decay_rows(ld_decay),
    )
    log_msg(f"done aggregating genome stats rep={args.rep}")

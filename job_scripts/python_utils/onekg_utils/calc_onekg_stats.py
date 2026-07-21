###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           calc_onekg_stats.py
###############################################################################

# overview: coordinate empirical folded-SFS, diversity, LD, and kinship stats.


##### set up ##################################################################
from collections import Counter, defaultdict
from pathlib import Path
import gzip
import math
from .read_onekg_sample_pops import read_onekg_sample_pops
from .read_tsv_rows import read_tsv_rows
from .scan_onekg_vcf import scan_onekg_vcf
from .write_stats_table import write_stats_table


##### internal functions #####################################################
'''
internal: sum folded one-dimensional SFS bins across chromosome tables.
'''
def _aggregate_folded_sfs_rows(rows):
    grouped = defaultdict(float)
    for row in rows:
        key = (
            int(row["rep"]),
            row["pop"],
            int(row["minor_allele_count"]),
        )
        grouped[key] += float(row["count"])
    return [
        {
            "rep": rep,
            "pop": pop,
            "minor_allele_count": minor_count,
            "count": count,
        }
        for (rep, pop, minor_count), count in sorted(grouped.items())
    ]

'''
internal: aggregate chromosome LD-decay summaries at genome scope.
'''
def _aggregate_ld_rows(rows):
    grouped = defaultdict(lambda: {"sum_r2": 0.0, "n_pairs": 0})
    for row in rows:
        key = (
            int(row["rep"]),
            row["pop"],
            int(row["distance_bin_bp"]),
        )
        grouped[key]["sum_r2"] += float(row["sum_r2"])
        grouped[key]["n_pairs"] += int(row["n_pairs"])
    output = []
    for (rep, pop, distance_bin), values in sorted(grouped.items()):
        num_pairs = values["n_pairs"]
        output.append(
            {
                "rep": rep,
                "pop": pop,
                "distance_bin_bp": distance_bin,
                "mean_r2": (
                    values["sum_r2"] / num_pairs
                    if num_pairs else math.nan
                ),
                "sum_r2": values["sum_r2"],
                "n_pairs": num_pairs,
            }
        )
    return output

'''
internal: aggregate chromosome pi and theta rows at genome scope by rebuilding
their numerators from values and family-specific spans.
'''
def _aggregate_pi_theta_rows(rows):
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

'''
internal: bin PLINK LD pairs into nonoverlapping genomic windows and distance
bins.
'''
def _bin_plink_ld_rows(
    pairs,
    rep,
    chrom,
    pop,
    window_size_bp=2_000_000,
    distance_bin_bp=5_000,
):
    grouped = defaultdict(lambda: {"sum_r2": 0.0, "n_pairs": 0})
    for pair in pairs:
        pos1_key = next(
            (key for key in ("POS_A", "POS1", "BP_A") if key in pair),
            None,
        )
        pos2_key = next(
            (key for key in ("POS_B", "POS2", "BP_B") if key in pair),
            None,
        )
        r2_key = next(
            (key for key in ("UNPHASED_R2", "R2") if key in pair),
            None,
        )
        if pos1_key is None or pos2_key is None or r2_key is None:
            raise ValueError("PLINK LD output lacks position or r2 columns")
        pos1 = int(pair[pos1_key])
        pos2 = int(pair[pos2_key])
        window_start = ((pos1 - 1) // window_size_bp) * window_size_bp
        if ((pos2 - 1) // window_size_bp) * window_size_bp != window_start:
            continue
        distance = pos2 - pos1
        if distance <= 0 or distance > window_size_bp:
            continue
        distance_bin = (
            ((distance - 1) // distance_bin_bp) + 1
        ) * distance_bin_bp
        key = (window_start, distance_bin)
        r2_value = float(pair[r2_key])
        grouped[key]["sum_r2"] += (
            r2_value if math.isfinite(r2_value) else 0.0
        )
        grouped[key]["n_pairs"] += 1
    output = []
    for (window_start, distance_bin), values in sorted(grouped.items()):
        output.append(
            {
                "rep": rep,
                "chrom": chrom,
                "pop": pop,
                "window_start": window_start,
                "window_end": window_start + window_size_bp,
                "distance_bin_bp": distance_bin,
                "mean_r2": values["sum_r2"] / values["n_pairs"],
                "sum_r2": values["sum_r2"],
                "n_pairs": values["n_pairs"],
            }
        )
    return output

'''
internal: build a folded 1D SFS from REF and ALT allele counts; here, we do not
use polarization/ancestral-state inference.
'''
def _build_folded_1d_sfs_rows(rep, chrom, pop, allele_counts):
    bins = Counter(min(ref_count, alt_count)
                   for ref_count, alt_count in allele_counts)
    return [
        {
            "rep": rep,
            "chrom": chrom,
            "pop": pop,
            "minor_allele_count": minor_count,
            "count": count,
        }
        for minor_count, count in sorted(bins.items())
    ]

'''
internal: build pi and Watterson theta rows using one explicit sequence denom.
'''
def _build_pi_theta_rows(
    rep,
    chrom,
    pop,
    allele_counts,
    span,
    mutation_rate,
    num_haplotypes,
):
    if span <= 0:
        raise ValueError("Sequence span must be positive")
    sample_sizes = {
        ref_count + alt_count for ref_count, alt_count in allele_counts
    }
    if sample_sizes and sample_sizes != {num_haplotypes}:
        raise ValueError("Allele counts must have one called sample size")
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

'''
internal: stream dictionary rows from a whitespace-delimited PLINK LD report.
Header names are normalized by removing PLINK's optional leading hash char.
'''
def _read_plink_ld_rows(path):
    with open(path, "r", encoding="utf-8") as in_file:
        header_line = next((line for line in in_file if line.strip()), None)
        if header_line is None:
            raise ValueError(f"Empty PLINK LD output {path}")
        headers = [header.lstrip("#") for header in header_line.split()]
        for line in in_file:
            if line.strip():
                yield dict(zip(headers, line.split()))


'''
internal: read the post-QC callable length for one chromosome from GRCh38.p14.
'''
def _read_onekg_chr_lengths(path, chrom):
    with open(path, "r", encoding="utf-8") as in_file:
        rows = [line.split() for line in in_file if line.strip()]
    if not rows:
        raise ValueError("Chromosome-length file is empty")
    headers = [value.lower() for value in rows[0]]
    chrom_names = ("chrom", "chr", "chromosome")
    qc_names = ("chr_len_after_qc", "length_after_qc", "callable_length")
    chrom_index = next(
        (headers.index(name) for name in chrom_names if name in headers),
        None,
    )
    qc_index = next(
        (headers.index(name) for name in qc_names if name in headers),
        None,
    )
    if None in (chrom_index, qc_index):
        raise ValueError(
            "Chromosome-length file needs chrom and chr_len_after_qc columns"
        )
    target = str(chrom).removeprefix("chr")
    for fields in rows[1:]:
        row_chrom = fields[chrom_index].removeprefix("chr")
        if row_chrom == target:
            span = float(fields[qc_index])
            if span <= 0:
                raise ValueError(
                    f"Chromosome {chrom} has invalid chr_len_after_qc"
                )
            return span
    raise ValueError(f"Chromosome {chrom} is absent from chromosome lengths")

'''
calculate the union span of BED intervals for one chromosome. Coordinates use
zero-based, half-open BED semantics, and adjacent intervals are merged.
'''
def _calc_bed_union_span(rows, chrom):
    target_chrom = str(chrom).lower().removeprefix("chr")
    intervals = []
    for line_number, raw_line in enumerate(rows, start=1):
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        fields = line.split()
        if len(fields) < 3:
            raise ValueError(
                f"Malformed BED interval on line {line_number}: "
                "expected at least three columns"
            )
        row_chrom = fields[0].lower().removeprefix("chr")
        try:
            start = int(fields[1])
            end = int(fields[2])
        except ValueError as error:
            raise ValueError(
                f"Malformed BED interval on line {line_number}: "
                "coordinates must be integers"
            ) from error
        if start < 0 or end < start:
            raise ValueError(
                f"Malformed BED interval on line {line_number}: "
                f"invalid half-open interval [{start}, {end})"
            )
        if row_chrom == target_chrom and end > start:
            intervals.append((start, end))

    if not intervals:
        raise ValueError(
            f"Chromosome {chrom} has no positive BED span"
        )

    intervals.sort()
    span = 0
    current_start, current_end = intervals[0]
    for start, end in intervals[1:]:
        if start <= current_end:
            current_end = max(current_end, end)
        else:
            span += current_end - current_start
            current_start, current_end = start, end
    return span + current_end - current_start

##### main function ###########################################################
'''
calculate chromosome statistics or aggregate genome statistics for one
population. Writes canonical diversity, folded-SFS, LD, KING, and QC tables.
'''
def calc_onekg_stats(args):
    stats_dir = Path(args.stats_dir)
    stats_dir.mkdir(parents=True, exist_ok=True)

    if args.analysis_level not in {"chromosome", "genome"}:
        raise ValueError(f"Unknown analysis level {args.analysis_level}")

    if args.analysis_level == "chromosome":
        required_args = {
            "ld_path": args.ld_path,
            "vcf_path": args.vcf_path,
            "intergenic_vcf_path": args.intergenic_vcf_path,
            "intergenic_bed_path": args.intergenic_bed_path,
            "unrels_path": args.unrels_path,
            "fam_path": args.fam_path,
            "chr_lens_path": args.chr_lens_path,
            "chrom": args.chrom,
            "mutation_rate": args.mutation_rate,
            "pops": args.pops,
        }
        missing_args = [
            name for name, value in required_args.items() if value is None
        ]
        if missing_args:
            missing_flags = ", ".join(
                f"--{name.replace('_', '-')}" for name in missing_args
            )
            raise ValueError(
                f"Missing chromosome arguments: {missing_flags}"
            )
    elif not args.chroms:
        raise ValueError("Missing genome argument: --chroms")

    # normalize KING output before writing the unrelated-kinship table.
    with open(args.king_path, "r", encoding="utf-8") as in_file:
        lines = [line.split() for line in in_file if line.strip()]
    if not lines:
        raise ValueError(f"Empty KING output {args.king_path}")
    headers = [header.lstrip("#") for header in lines[0]]
    king_rows = []
    for values in lines[1:]:
        raw = dict(zip(headers, values))
        row = {"rep": 0}
        if args.analysis_level == "chromosome":
            row["chrom"] = args.chrom
        row.update(
            {
                "pop": args.pop,
                "id1": raw.get("IID1", raw.get("ID1")),
                "id2": raw.get("IID2", raw.get("ID2")),
                "kinship": raw.get("KINSHIP", raw.get("Kinship")),
            }
        )
        king_rows.append(row)
    king_suffix = (
        f"rep_0.chr{args.chrom}.{args.pop}"
        if args.analysis_level == "chromosome"
        else f"rep_0.{args.pop}"
    )
    write_stats_table(
        stats_dir / f"kinship_unrelated.{king_suffix}",
        king_rows,
    )

    # chromosome jobs scan genotypes and calculate population statistics.
    if args.analysis_level == "chromosome":
        sample_pops = read_onekg_sample_pops(
            args.unrels_path,
            args.fam_path,
            args.pops,
        )
        open_vcf = gzip.open if args.vcf_path.endswith(".gz") else open
        with open_vcf(args.vcf_path, "rt", encoding="utf-8") as vcf_file:
            scan = scan_onekg_vcf(vcf_file, sample_pops, args.pops)

        pop_counts = scan["counts_by_pop"][args.pop]
        allele_counts = list(pop_counts.values())
        num_haplotypes = 2 * scan["sample_counts"][args.pop]
        callable_span = _read_onekg_chr_lengths(
            args.chr_lens_path,
            args.chrom,
        )
        intergenic_sample_pops = {
            sample: pop
            for sample, pop in sample_pops.items()
            if pop == args.pop
        }
        open_intergenic_vcf = (
            gzip.open
            if args.intergenic_vcf_path.endswith(".gz")
            else open
        )
        with open_intergenic_vcf(
            args.intergenic_vcf_path,
            "rt",
            encoding="utf-8",
        ) as vcf_file:
            intergenic_scan = scan_onekg_vcf(
                vcf_file,
                intergenic_sample_pops,
                [args.pop],
            )
        intergenic_counts = list(
            intergenic_scan["counts_by_pop"][args.pop].values()
        )
        with open(
            args.intergenic_bed_path,
            "r",
            encoding="utf-8",
        ) as bed_file:
            intergenic_span = _calc_bed_union_span(
                bed_file,
                args.chrom,
            )
        pi_theta_intergenic = _build_pi_theta_rows(
            0,
            args.chrom,
            args.pop,
            intergenic_counts,
            intergenic_span,
            args.mutation_rate,
            num_haplotypes,
        )
        pi_theta_callable = _build_pi_theta_rows(
            0,
            args.chrom,
            args.pop,
            allele_counts,
            callable_span,
            args.mutation_rate,
            num_haplotypes,
        )
        sfs = _build_folded_1d_sfs_rows(
            0,
            args.chrom,
            args.pop,
            allele_counts,
        )
        counts_rows = [
            {
                "rep": 0,
                "chrom": site.split(":", 1)[0],
                "position": int(site.split(":", 1)[1]),
                "pop": args.pop,
                "ref_count": counts[0],
                "alt_count": counts[1],
            }
            for site, counts in sorted(
                pop_counts.items(),
                key=lambda item: int(item[0].split(":", 1)[1]),
            )
        ]
        qc_row = {
            "rep": 0,
            "chrom": args.chrom,
            "pop": args.pop,
            **scan["qc"],
            **{
                f"retained_{pop}_samples": scan["sample_counts"][pop]
                for pop in args.pops
            },
        }
        outputs = {
            (
                f"allele_counts.rep_0.chr{args.chrom}.{args.pop}"
            ): counts_rows,
            (
                f"pi_theta_stats_intergenic.rep_0.chr{args.chrom}."
                f"{args.pop}"
            ): pi_theta_intergenic,
            (
                f"pi_theta_stats_full_callable_chrom.rep_0."
                f"chr{args.chrom}."
                f"{args.pop}"
            ): pi_theta_callable,
            f"sfs.rep_0.chr{args.chrom}.{args.pop}": sfs,
            (
                f"variant_qc.rep_0.chr{args.chrom}.{args.pop}"
            ): [qc_row],
            (
                f"ld_decay.rep_0.chr{args.chrom}.{args.pop}"
            ): _bin_plink_ld_rows(
                _read_plink_ld_rows(args.ld_path),
                0,
                args.chrom,
                args.pop,
            ),
        }
        for output_name, rows in outputs.items():
            write_stats_table(stats_dir / output_name, rows)
        return

    for table_name in (
        "pi_theta_stats_intergenic",
        "pi_theta_stats_full_callable_chrom",
    ):
        rows = []
        for chrom in args.chroms:
            path = stats_dir / (
                f"{table_name}.rep_0.chr{chrom}.{args.pop}.tsv"
            )
            rows.extend(read_tsv_rows(path))
        aggregated = _aggregate_pi_theta_rows(rows)
        write_stats_table(
            stats_dir / f"{table_name}.rep_0.{args.pop}",
            aggregated,
        )

    sfs_rows = []
    ld_rows = []
    for chrom in args.chroms:
        sfs_rows.extend(
            read_tsv_rows(
                stats_dir / f"sfs.rep_0.chr{chrom}.{args.pop}.tsv"
            )
        )
        ld_rows.extend(
            read_tsv_rows(
                stats_dir / f"ld_decay.rep_0.chr{chrom}.{args.pop}.tsv"
            )
        )
    write_stats_table(
        stats_dir / f"sfs.rep_0.{args.pop}",
        _aggregate_folded_sfs_rows(sfs_rows),
    )
    write_stats_table(
        stats_dir / f"ld_decay.rep_0.{args.pop}",
        _aggregate_ld_rows(ld_rows),
    )

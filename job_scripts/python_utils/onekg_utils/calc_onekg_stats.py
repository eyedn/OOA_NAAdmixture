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
from pathlib import Path
import gzip
from .aggregate_folded_sfs_rows import aggregate_folded_sfs_rows
from .aggregate_ld_rows import aggregate_ld_rows
from .aggregate_pi_theta_rows import aggregate_pi_theta_rows
from .bin_plink_ld_rows import bin_plink_ld_rows
from .build_folded_1d_sfs_rows import build_folded_1d_sfs_rows
from .build_pi_theta_rows import build_pi_theta_rows
from .read_plink_ld_rows import read_plink_ld_rows
from .read_onekg_chr_lengths import read_onekg_chr_lengths
from .read_onekg_sample_pops import read_onekg_sample_pops
from .read_tsv_rows import read_tsv_rows
from .scan_onekg_vcf import scan_onekg_vcf
from .write_stats_table import write_stats_table


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
        chr_len_after_qc, chr_len = read_onekg_chr_lengths(
            args.chr_lens_path,
            args.chrom,
        )
        pi_theta = build_pi_theta_rows(
            0,
            args.chrom,
            args.pop,
            allele_counts,
            chr_len_after_qc,
            args.mutation_rate,
        )
        pi_theta_full = build_pi_theta_rows(
            0,
            args.chrom,
            args.pop,
            allele_counts,
            chr_len,
            args.mutation_rate,
        )
        sfs = build_folded_1d_sfs_rows(
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
                f"pi_theta_stats.rep_0.chr{args.chrom}.{args.pop}"
            ): pi_theta,
            (
                f"pi_theta_stats_full_chrom.rep_0.chr{args.chrom}."
                f"{args.pop}"
            ): pi_theta_full,
            f"sfs.rep_0.chr{args.chrom}.{args.pop}": sfs,
            (
                f"variant_qc.rep_0.chr{args.chrom}.{args.pop}"
            ): [qc_row],
            (
                f"ld_decay.rep_0.chr{args.chrom}.{args.pop}"
            ): bin_plink_ld_rows(
                read_plink_ld_rows(args.ld_path),
                0,
                args.chrom,
                args.pop,
            ),
        }
        for output_name, rows in outputs.items():
            write_stats_table(stats_dir / output_name, rows)
        return

    haplotype_counts = set()
    for chrom in args.chroms:
        counts_path = stats_dir / (
            f"allele_counts.rep_0.chr{chrom}.{args.pop}.tsv"
        )
        counts_rows = read_tsv_rows(counts_path)
        if not counts_rows:
            raise ValueError(f"No allele counts in {counts_path}")
        haplotype_counts.add(
            int(counts_rows[0]["ref_count"])
            + int(counts_rows[0]["alt_count"])
        )
    if len(haplotype_counts) != 1:
        raise ValueError(
            f"Inconsistent haplotype counts for pop={args.pop}: "
            f"{sorted(haplotype_counts)}"
        )
    haplotypes_by_pop = {args.pop: haplotype_counts.pop()}

    for table_name in (
        "pi_theta_stats",
        "pi_theta_stats_full_chrom",
    ):
        rows = []
        for chrom in args.chroms:
            path = stats_dir / (
                f"{table_name}.rep_0.chr{chrom}.{args.pop}.tsv"
            )
            rows.extend(read_tsv_rows(path))
        aggregated = aggregate_pi_theta_rows(rows, haplotypes_by_pop)
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
        aggregate_folded_sfs_rows(sfs_rows),
    )
    write_stats_table(
        stats_dir / f"ld_decay.rep_0.{args.pop}",
        aggregate_ld_rows(ld_rows),
    )

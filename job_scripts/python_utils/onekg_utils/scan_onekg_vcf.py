###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           scan_onekg_vcf.py
###############################################################################

# overview: scan VCF records into population-aware genotype and QC summaries.


##### set up ##################################################################
from .parse_gt import parse_gt


##### main function ###########################################################
'''
scan a VCF and retain complete biallelic sites across requested samples.
Returns population allele counts, retained sample counts, and QC totals.
'''
def scan_onekg_vcf(vcf_file, sample_pops, pops):
    sample_names = None
    sample_indexes = None
    counts_by_pop = {pop: {} for pop in pops}
    sample_counts = {
        pop: sum(sample_pop == pop for sample_pop in sample_pops.values())
        for pop in pops
    }
    qc = {
        "input_variant_count": 0,
        "non_biallelic_or_malformed_variant_count": 0,
        "missing_genotype_removal_count": 0,
        "complete_site_retained_count": 0,
        "folded_sfs_variant_count": 0,
    }

    for raw_line in vcf_file:
        # skip metadata until the sample-bearing column header is reached.
        if raw_line.startswith("##"):
            continue
        if raw_line.startswith("#CHROM"):
            sample_names = raw_line.rstrip("\n").split("\t")[9:]
            missing = sorted(set(sample_pops) - set(sample_names))
            if missing:
                raise ValueError(
                    "Requested samples missing from VCF: " + ", ".join(missing)
                )
            sample_indexes = [
                (index, sample_pops[sample])
                for index, sample in enumerate(sample_names)
                if sample in sample_pops
            ]
            continue
        if raw_line.startswith("#") or not raw_line.strip():
            continue
        if sample_indexes is None:
            raise ValueError("VCF header is missing #CHROM")

        # reject malformed, multiallelic, or incompletely called sites.
        qc["input_variant_count"] += 1
        fields = raw_line.rstrip("\n").split("\t")
        if len(fields) < 10 or "," in fields[4] or fields[4] == ".":
            qc["non_biallelic_or_malformed_variant_count"] += 1
            continue
        format_fields = fields[8].split(":")
        if "GT" not in format_fields:
            qc["non_biallelic_or_malformed_variant_count"] += 1
            continue
        gt_index = format_fields.index("GT")
        site_counts = {pop: [0, 0] for pop in pops}
        missing_gt = False
        malformed_gt = False
        for sample_index, pop in sample_indexes:
            sample_fields = fields[9 + sample_index].split(":")
            if gt_index >= len(sample_fields):
                malformed_gt = True
                break
            try:
                counts = parse_gt(sample_fields[gt_index])
            except ValueError:
                malformed_gt = True
                break
            if counts is None:
                missing_gt = True
                break
            site_counts[pop][0] += counts[0]
            site_counts[pop][1] += counts[1]
        if malformed_gt:
            qc["non_biallelic_or_malformed_variant_count"] += 1
            continue
        if missing_gt:
            qc["missing_genotype_removal_count"] += 1
            continue

        # retain polymorphic complete sites for folded frequency summaries.
        qc["complete_site_retained_count"] += 1
        total_ref = sum(counts[0] for counts in site_counts.values())
        total_alt = sum(counts[1] for counts in site_counts.values())
        if total_ref > 0 and total_alt > 0:
            site_key = f"{fields[0]}:{fields[1]}"
            for pop in pops:
                counts_by_pop[pop][site_key] = tuple(site_counts[pop])
            qc["folded_sfs_variant_count"] += 1

    if sample_names is None:
        raise ValueError("VCF header is missing #CHROM")
    return {
        "counts_by_pop": counts_by_pop,
        "sample_counts": sample_counts,
        "qc": qc,
    }

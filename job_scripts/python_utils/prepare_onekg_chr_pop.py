###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           prepare_onekg_chr_pop.py
###############################################################################

# overview: prepare empirical chrom.-population subset for downstream tools.


##### set up ##################################################################
from pathlib import Path
import argparse
import gzip
from onekg_utils.read_onekg_sample_pops import read_onekg_sample_pops
from onekg_utils.scan_onekg_vcf import scan_onekg_vcf


##### arguments ###############################################################
'''
define the source VCF and metadata paths, output directory, chromosome,
target population, and ordered population list.
'''
parser = argparse.ArgumentParser()
parser.add_argument("--vcf-path", required=True)
parser.add_argument("--unrels-path", required=True)
parser.add_argument("--fam-path", required=True)
parser.add_argument("--output-dir", required=True)
parser.add_argument("--chrom", required=True)
parser.add_argument("--pop", required=True)
parser.add_argument("--pops", nargs="+", required=True)


##### main ####################################################################
if __name__ == "__main__":
    args = parser.parse_args()
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    sample_pops = read_onekg_sample_pops(
        args.unrels_path,
        args.fam_path,
        args.pops,
    )
    open_vcf = gzip.open if args.vcf_path.endswith(".gz") else open

    # recover source VCF order so every sample handoff remains deterministic.
    with open_vcf(args.vcf_path, "rt", encoding="utf-8") as vcf_file:
        scan = scan_onekg_vcf(vcf_file, sample_pops, args.pops)

    with open_vcf(args.vcf_path, "rt", encoding="utf-8") as vcf_file:
        for line in vcf_file:
            if line.startswith("#CHROM"):
                vcf_samples = line.rstrip("\n").split("\t")[9:]
                break

    prefix = f"onekg.rep_0.chr{args.chrom}"
    ordered_samples = [
        sample for sample in vcf_samples if sample in sample_pops
    ]

    # these path should matches their analogous variables in the bash script
    all_samples_path = output_dir / f"{prefix}.all.{args.pop}.samples.tmp"
    pop_samples_path = output_dir / f"{prefix}.{args.pop}.samples"
    pop_keep_path = output_dir / f"{prefix}.{args.pop}.keep"
    sites_path = output_dir / f"{prefix}.complete_sites.{args.pop}.tsv"

    # write all-sample, population, PLINK keep, and complete-site handoffs.
    with open(all_samples_path, "w", encoding="utf-8") as out_file:
        out_file.writelines(f"{sample}\n" for sample in ordered_samples)
    pop_samples = [
        sample for sample in ordered_samples if sample_pops[sample] == args.pop
    ]
    if not pop_samples:
        raise ValueError(f"No retained samples found for pop={args.pop}")
    with open(pop_samples_path, "w", encoding="utf-8") as out_file:
        out_file.writelines(f"{sample}\n" for sample in pop_samples)
    # note, the format f"0\t{sample}\n" is due to a vcf being used as the 
    # starting position, where do FID is provided 
    with open(pop_keep_path, "w", encoding="utf-8") as out_file:
        out_file.writelines(f"0\t{sample}\n" for sample in pop_samples)
    retained_sites = scan["counts_by_pop"][args.pops[0]]
    with open(sites_path, "w", encoding="utf-8") as out_file:
        for site in retained_sites:
            chrom, position = site.split(":", 1)
            out_file.write(f"{chrom}\t{position}\n")

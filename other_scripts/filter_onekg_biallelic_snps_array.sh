#!/usr/bin/env bash

###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           filter_onekg_biallelic_snps_array.sh
###############################################################################

#SBATCH --cpus-per-task=4
#SBATCH --array=1-22
#SBATCH --mem=16
#SBATCH --time=1-00:01:00
#SBATCH --partition=qcb
#SBATCH --account=jazlynmo_738
#SBATCH --job-name="filterUnmasked2kG"
#SBATCH --nodes=1
#SBATCH --output="/home1/karatas/logs/tmp/%A_%a.%x.out"
#SBATCH --error="/home1/karatas/logs/tmp/%A_%a.%x.err"
#SBATCH --mail-type=ALL
#SBATCH --mail-user=karatas@usc.edu

# create one indexed, biallelic-SNP 1kG VCF for each chromosome array task.


##### set up ##################################################################
set -euo pipefail

: "${SLURM_ARRAY_TASK_ID:?ERROR: run with sbatch --array=1-22}"

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
project_dir="/home1/karatas/scratch/OOA_NAAdmixture_1kGwoStrickMask"
chrom="${SLURM_ARRAY_TASK_ID}"

module load htslib/1.19.1 bcftools/1.19


##### input and output paths ##################################################
samples_keep="${project_dir}/samples.keep"
source_dir="${HOME}/1000GenomeNYGC_hg38/vcfs"
source_stem="CCDG_14151_B01_GRM_WGS_2020-08-05_chr${chrom}"
source_suffix=".filtered.shapeit2-duohmm-phased.nodupmarkers.snps.vcf.gz"
source_vcf="${source_dir}/${source_stem}${source_suffix}"
output_dir="${project_dir}/vcfs"
output_vcf="${output_dir}/chr${chrom}.vcf.gz"

if [[ ! -s "${samples_keep}" ]]; then
    echo "ERROR: missing or empty sample keep file: ${samples_keep}" >&2
    exit 1
fi

if [[ ! -s "${source_vcf}" ]]; then
    echo "ERROR: missing or empty source VCF: ${source_vcf}" >&2
    exit 1
fi


##### biallelic-SNP filtering #################################################
mkdir -p "${output_dir}"

# subset first so the minor allele count reflects only retained samples.
bcftools view \
    --threads 4 \
    --samples-file "${samples_keep}" \
    --output-type u \
    "${source_vcf}" | \
bcftools view \
    --threads 4 \
    --min-alleles 2 \
    --max-alleles 2 \
    --types snps \
    --min-ac 1:minor \
    --output-type z \
    --output-file "${output_vcf}"
tabix -f -p vcf "${output_vcf}"

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

# create one indexed, biallelic-SNP 1kG VCF for each chromosome array task.


##### set up ##################################################################
set -euo pipefail

: "${SLURM_ARRAY_TASK_ID:?ERROR: run with sbatch --array=1-22}"

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
project_dir="$(cd "${script_dir}/.." && pwd)"
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

bcftools view \
    --threads 4 \
    --samples-file "${samples_keep}" \
    --min-alleles 2 \
    --max-alleles 2 \
    --types snps \
    --output-type z \
    --output-file "${output_vcf}" \
    "${source_vcf}"
tabix -f -p vcf "${output_vcf}"

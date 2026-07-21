#!/usr/bin/env bash

###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           parse_onekg_chr_pop.sh
###############################################################################

# workflow: prepare one 1000 Genomes chromosome-population analysis input.


##### set up ##################################################################
set -euo pipefail

# load the software environment used for VCF, PLINK, and Python preparation.
module purge
ml gcc/13.3.0 htslib/1.19.1 bcftools/1.19 plink2/2.00a4.3 conda
source /apps/conda/miniforge3/25.3.0/etc/profile.d/conda.sh
conda activate OOA_NAAdmixture

# load shared logging and chromosome-population task mapping.
project_dir="$(pwd)"
source "${project_dir}/other_scripts/log_msg.sh"
source "${project_dir}/other_scripts/map_slurm_task.sh"
: "${SLURM_ARRAY_TASK_ID:?ERROR: run as a Slurm array}"


##### variables ###############################################################
# read fixed inputs followed by chromosome and population arrays.
vcf_prefix="$1"
vcf_suffix="$2"
unrels_path="$3"
fam_path="$4"
out_vcf_dir="$5"
out_bed_dir="$6"
pop_info_dir="$7"
shift 7
shift   # skip the "--" from input arguments
chroms=()
while [[ "$1" != "--" ]]; do
    chroms+=( "$1" )
    shift
done
shift   # skip the "--" from input arguments
pops=( "$@" )

# map this chromosome-major task ID and derive all output prefixes.
read -r chr pop < <(
    map_slurm_task_to_chr_pop \
        "${SLURM_ARRAY_TASK_ID}" \
        -- \
        "${chroms[@]}" \
        -- \
        "${pops[@]}"
)
num_threads="${SLURM_CPUS_PER_TASK:-1}"
input_vcf="${vcf_prefix}${chr}${vcf_suffix}"
prefix="onekg.rep_0.chr${chr}"
common_vcf="${out_vcf_dir}/${prefix}.all.${pop}.vcf.gz"
pop_vcf="${out_vcf_dir}/${prefix}.${pop}.vcf.gz"
pop_bed_prefix="${out_bed_dir}/${prefix}.${pop}"


##### sample and variant filtering ############################################
# create output directories and population-specific sample and site lists.
mkdir -p "${out_vcf_dir}" "${out_bed_dir}" "${pop_info_dir}"

# write population sample lists, a PLINK keep file, and complete shared sites.
python "${project_dir}/job_scripts/python_utils/prepare_onekg_chr_pop.py" \
    --vcf-path "${input_vcf}" \
    --unrels-path "${unrels_path}" \
    --fam-path "${fam_path}" \
    --output-dir "${pop_info_dir}" \
    --chrom "${chr}" \
    --pop "${pop}" \
    --pops "${pops[@]}"
# these path should matches their analogous variables in the above py. script
all_samples="${pop_info_dir}/${prefix}.all.${pop}.samples.tmp"
pop_samples="${pop_info_dir}/${prefix}.${pop}.samples"
complete_sites="${pop_info_dir}/${prefix}.complete_sites.${pop}.tsv"

# filter the source VCF to complete shared sites and the selected population.
log_msg "filtering complete biallelic sites chr=${chr} pop=${pop}"
bcftools view \
    --threads "${num_threads}" \
    -S "${all_samples}" \
    -R "${complete_sites}" \
    -m2 -M2 -v snps -g ^miss \
    -Oz -o "${common_vcf}" \
    "${input_vcf}"
tabix -f -p vcf "${common_vcf}"
rm -rf "$all_samples"

bcftools view \
    --threads "${num_threads}" \
    -S "${pop_samples}" \
    -Oz -o "${pop_vcf}" \
    "${common_vcf}"
tabix -f -p vcf "${pop_vcf}"


##### PLINK output ############################################################
# generate PLINK binaries for the population-specific VCF.
plink2 \
    --vcf "${pop_vcf}" \
    --set-all-var-ids '@:#:$r:$a' \
    --threads "${num_threads}" \
    --make-bed \
    --out "${pop_bed_prefix}"

# retain one chromosome-level common-site reference without duplicate writers.
if [[ "${pop}" == "${pops[0]}" ]]; then
    cp "${complete_sites}" \
        "${pop_info_dir}/common_biallelic_complete.chr${chr}.tsv"
fi
log_msg "done parsing 1000 Genomes chr=${chr} pop=${pop}"

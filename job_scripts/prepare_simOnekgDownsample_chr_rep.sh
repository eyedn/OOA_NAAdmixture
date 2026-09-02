#!/usr/bin/env bash

###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           prepare_simOnekgDownsample_chr_rep.sh
###############################################################################

# workflow: match one simulation chromosome to empirical 10-kb SNP density.


##### set up ##################################################################
set -euo pipefail

conda_env="$1"
shift

module purge
ml gcc/13.3.0 htslib/1.19.1 bcftools/1.19 vcftools/0.1.16 \
    plink2/2.00a4.3 conda
source /apps/conda/miniforge3/25.3.0/etc/profile.d/conda.sh
conda activate "${conda_env}"

project_dir="$(pwd)"
source "${project_dir}/other_scripts/log_msg.sh"
source "${project_dir}/other_scripts/map_slurm_task.sh"
: "${SLURM_ARRAY_TASK_ID:?ERROR: run as a Slurm array}"


##### variables ###############################################################
# read fixed inputs followed by chromosome and population arrays.
sim_vcf_dir="$1"
sim_pop_info_dir="$2"
empirical_vcf_dir="$3"
empirical_admixed_pop="$4"
span_incl_file="$5"
intergenic_file="$6"
simOnekgDownsample_vcf_dir="$7"
simOnekgDownsample_plink_bed_dir="$8"
simOnekgDownsample_pop_info_dir="$9"
simOnekgDownsample_admixture_dir="${10}"
simOnekgDownsample_faststructure_dir="${11}"
simOnekgDownsample_king_dir="${12}"
simOnekgDownsample_stats_dir="${13}"
snp_density_window_bp="${14}"
num_reps="${15}"
genetic_map="${16}"
shift 16
shift 1 # skip the "--" from input arguments
chroms=()
while [[ "$1" != "--" ]]; do
    chroms+=( "$1" )
    shift
done
shift 1 # skip the "--" from input arguments
pops=( "$@" )

read -r chr rep < <(
    map_slurm_task_to_chr_rep \
        "${SLURM_ARRAY_TASK_ID}" \
        "${num_reps}" \
        -- \
        "${chroms[@]}"
)

# derive task-specific input, output, and intermediate paths.
num_threads="${SLURM_CPUS_PER_TASK:-1}"
prefix="${genetic_map}_${rep}_chr${chr}_all"
metadata_prefix="${genetic_map}_${rep}_chr${chr}"
sim_vcf="${sim_vcf_dir}/${prefix}.biallelic_snps.vcf.gz"
callable_sim_vcf="${simOnekgDownsample_vcf_dir}/${prefix}.strict_callable.vcf.gz"
empirical_vcf="${empirical_vcf_dir}/onekg.rep_0.chr${chr}.${empirical_admixed_pop}.vcf.gz"
simOnekgDownsample_vcf="${simOnekgDownsample_vcf_dir}/${prefix}.biallelic_snps.vcf.gz"
sim_metadata="${sim_pop_info_dir}/${metadata_prefix}.sample_metadata.tsv"
simOnekgDownsample_metadata="${simOnekgDownsample_pop_info_dir}/${metadata_prefix}.sample_metadata.tsv"
audit_prefix="${simOnekgDownsample_pop_info_dir}/${prefix}.snp_density"
empirical_density_prefix="${audit_prefix}.empirical"
callable_density_prefix="${audit_prefix}.simulation_callable"
sim_positions="${simOnekgDownsample_pop_info_dir}/${prefix}.simulation_positions.tsv"
wanted_snps="${simOnekgDownsample_pop_info_dir}/${prefix}.wanted_snps.tsv"
contract_json="${simOnekgDownsample_pop_info_dir}/${prefix}.snp_density.contract.json"
selection_seed=$((1000 * rep + chr))


##### density downsampling ####################################################
mkdir -p \
    "${simOnekgDownsample_vcf_dir}" \
    "${simOnekgDownsample_plink_bed_dir}" \
    "${simOnekgDownsample_pop_info_dir}" \
    "${simOnekgDownsample_admixture_dir}" \
    "${simOnekgDownsample_faststructure_dir}" \
    "${simOnekgDownsample_king_dir}" \
    "${simOnekgDownsample_stats_dir}"

for required_path in \
    "${sim_vcf}" \
    "${empirical_vcf}" \
    "${sim_metadata}" \
    "${span_incl_file}" \
    "${intergenic_file}"; do
    if [[ ! -s "${required_path}" ]]; then
        echo "ERROR: missing required input ${required_path}" >&2
        exit 1
    fi
done

log_msg "calculating SNP density rep=${rep} chr=${chr}"
bcftools view \
    --threads "${num_threads}" \
    --regions-file "${span_incl_file}" \
    --output-type z \
    --output "${callable_sim_vcf}" \
    "${sim_vcf}"
tabix -f -p vcf "${callable_sim_vcf}"
vcftools \
    --gzvcf "${empirical_vcf}" \
    --SNPdensity "${snp_density_window_bp}" \
    --out "${empirical_density_prefix}"
vcftools \
    --gzvcf "${callable_sim_vcf}" \
    --SNPdensity "${snp_density_window_bp}" \
    --out "${callable_density_prefix}"
bcftools query \
    --format '%CHROM\t%POS\n' \
    "${callable_sim_vcf}" > "${sim_positions}"

python \
"${project_dir}"/job_scripts/python_utils/select_simOnekgDownsample_snps.py \
    --empirical-density "${empirical_density_prefix}.snpden" \
    --simulation-density "${callable_density_prefix}.snpden" \
    --simulation-positions "${sim_positions}" \
    --wanted-snps "${wanted_snps}" \
    --contract-json "${contract_json}" \
    --stats-dir "${simOnekgDownsample_stats_dir}" \
    --window-size-bp "${snp_density_window_bp}" \
    --seed "${selection_seed}" \
    --rep "${rep}" \
    --chrom "${chr}"
if [[ ! -s "${wanted_snps}" ]]; then
    echo "ERROR: SNP-density selection retained no variants" >&2
    exit 1
fi

bcftools view \
    --threads "${num_threads}" \
    --regions-file "${wanted_snps}" \
    --output-type z \
    --output "${simOnekgDownsample_vcf}" \
    "${callable_sim_vcf}"
tabix -f -p vcf "${simOnekgDownsample_vcf}"
if [[ ! -s "${simOnekgDownsample_vcf}" \
    || ! -s "${simOnekgDownsample_vcf}.tbi" ]]; then
    echo "ERROR: failed to create indexed VCF ${simOnekgDownsample_vcf}" >&2
    exit 1
fi


##### population inputs ######################################################
cp "${sim_metadata}" "${simOnekgDownsample_metadata}"
plink2 \
    --vcf "${simOnekgDownsample_vcf}" \
    --set-all-var-ids '@:#:$r:$a' \
    --threads "${num_threads}" \
    --make-bed \
    --out "${simOnekgDownsample_plink_bed_dir}/${prefix}"
for extension in bed bim fam; do
    bed_path="${simOnekgDownsample_plink_bed_dir}/${prefix}.${extension}"
    if [[ ! -s "${bed_path}" ]]; then
        echo "ERROR: missing all-population PLINK file ${bed_path}" >&2
        exit 1
    fi
done

for pop in "${pops[@]}"; do
    pop_prefix="${metadata_prefix}_${pop}"
    sample_path="${simOnekgDownsample_pop_info_dir}/${pop_prefix}.samples"
    keep_path="${simOnekgDownsample_pop_info_dir}/${pop_prefix}.keep"
    pop_vcf="${simOnekgDownsample_vcf_dir}/${pop_prefix}.vcf.gz"
    intergenic_vcf="${simOnekgDownsample_vcf_dir}/${pop_prefix}.intergenic.vcf.gz"
    pop_bed_prefix="${simOnekgDownsample_plink_bed_dir}/${pop_prefix}"

    awk -F '\t' -v pop="${pop}" 'NR > 1 && $3 == pop {print $2}' \
        "${simOnekgDownsample_metadata}" > "${sample_path}"
    awk -F '\t' -v pop="${pop}" \
        'NR > 1 && $3 == pop {print $1 "\t" $2}' \
        "${simOnekgDownsample_metadata}" > "${keep_path}"
    if [[ ! -s "${sample_path}" || ! -s "${keep_path}" ]]; then
        echo "ERROR: no metadata samples found for pop=${pop}" >&2
        exit 1
    fi

    bcftools view \
        --threads "${num_threads}" \
        --samples-file "${sample_path}" \
        --output-type z \
        --output "${pop_vcf}" \
        "${simOnekgDownsample_vcf}"
    tabix -f -p vcf "${pop_vcf}"
    bcftools view \
        --threads "${num_threads}" \
        --regions-file "${intergenic_file}" \
        --output-type z \
        --output "${intergenic_vcf}" \
        "${pop_vcf}"
    tabix -f -p vcf "${intergenic_vcf}"
    plink2 \
        --vcf "${pop_vcf}" \
        --set-all-var-ids '@:#:$r:$a' \
        --threads "${num_threads}" \
        --make-bed \
        --out "${pop_bed_prefix}"
    for extension in bed bim fam; do
        if [[ ! -s "${pop_bed_prefix}.${extension}" ]]; then
            echo "ERROR: missing ${pop} PLINK ${extension} output" >&2
            exit 1
        fi
    done
done

log_msg "done preparing simOnekgDownsample rep=${rep} chr=${chr}"

#!/usr/bin/env bash

###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           sim_select_stdpopsim_models_worker.sh
###############################################################################

# workflow: simulate and export one selected stdpopsim model on chrom. 22.


##### set up ##################################################################
set -euo pipefail

: "${SLURM_ARRAY_TASK_ID:?ERROR: run as a Slurm array task}"
if (( $# != 4 )); then
    echo "ERROR: expected conda env, output dir, map, and msprime model" >&2
    exit 1
fi

conda_env="$1"
select_out_dir="$2"
genetic_map="$3"
msprime_model="$4"

module purge
ml gcc/13.3.0 htslib/1.19.1 bcftools/1.19 conda
source /apps/conda/miniforge3/25.3.0/etc/profile.d/conda.sh
conda activate "${conda_env}"
export PATH="${HOME}/.conda/envs/${conda_env}/bin:${PATH}"

project_dir="$(pwd)"
source "${project_dir}/other_scripts/log_msg.sh"


##### model samples ###########################################################
case "${SLURM_ARRAY_TASK_ID}" in
    1)
        model="OutOfAfrica_3G09"
        samples=(YRI:500 CEU:500 CHB:500)
        expected_samples=1500
        ;;
    2)
        model="OutOfAfrica_2T12"
        samples=(AFR:500 EUR:500)
        expected_samples=1000
        ;;
    3)
        model="AmericanAdmixture_4B18"
        samples=(AFR:500 EUR:500 ASIA:500 ADMIX:500)
        expected_samples=2000
        ;;
    *)
        echo "ERROR: invalid array task ID: ${SLURM_ARRAY_TASK_ID}" >&2
        exit 1
        ;;
esac

seed=1234
prefix="${genetic_map}_1_chr22_${model}"
tree_dir="${select_out_dir}/${model}/trees"
vcf_dir="${select_out_dir}/${model}/vcfs"
tree_path="${tree_dir}/${prefix}.ts"
compressed_tree_path="${tree_path}.tsz"
vcf_path="${vcf_dir}/${prefix}.vcf.gz"
mkdir -p "${tree_dir}" "${vcf_dir}"


##### simulation ##############################################################
log_msg "simulating ${model} chr22 with seed=${seed}"
stdpopsim -e msprime --msprime-model "${msprime_model}" \
    HomSap -c chr22 -g "${genetic_map}" -s "${seed}" \
    -d "${model}" -o "${tree_path}" "${samples[@]}"
if [[ ! -s "${tree_path}" ]]; then
    echo "ERROR: missing or empty tree sequence: ${tree_path}" >&2
    exit 1
fi


##### vcf and tree sequence ###################################################
log_msg "exporting and indexing VCF for ${model}"
python -m tskit vcf --contig-id chr22 "${tree_path}" | \
    bgzip -c > "${vcf_path}"
tabix -f -p vcf "${vcf_path}"

actual_samples=$(bcftools query -l "${vcf_path}" | wc -l)
if (( actual_samples != expected_samples )); then
    echo "ERROR: ${model} has ${actual_samples} VCF samples; " \
        "expected ${expected_samples}" >&2
    exit 1
fi
if [[ ! -s "${vcf_path}" || ! -s "${vcf_path}.tbi" ]]; then
    echo "ERROR: missing or empty indexed VCF: ${vcf_path}" >&2
    exit 1
fi

log_msg "compressing and verifying tree sequence for ${model}"
tszip -f -k "${tree_path}"
if [[ ! -s "${compressed_tree_path}" ]]; then
    echo "ERROR: missing compressed tree sequence: ${compressed_tree_path}" >&2
    exit 1
fi
tsunzip -c "${compressed_tree_path}" > /dev/null
rm "${tree_path}"
log_msg "completed ${model}; VCF samples=${actual_samples}"

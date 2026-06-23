#!/usr/bin/env bash

###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology 
#           Mooney Lab
#           ---
#           sim_model.sh
###############################################################################

#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --partition=qcb
#SBATCH --account=jazlynmo_738
#SBATCH --nodes=1
#SBATCH --output=/home1/karatas/logs/sim_OOANAA.%A_%a.%x.out
#SBATCH --error=/home1/karatas/logs/sim_OOANAA.%A_%a.%x.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=karatas@usc.edu


set -euo pipefail

module purge
ml gcc/13.3.0 htslib/1.19.1 plink2/2.00a4.3 conda
source /apps/conda/miniforge3/25.3.0/etc/profile.d/conda.sh
conda activate OOA_NAAdmixture
export PATH="${HOME}/.conda/envs/OOA_NAAdmixture/bin:${PATH}"

project_dir="$(pwd)"
source "${project_dir}/other_scripts/log_msg.sh"

: "${SLURM_ARRAY_TASK_ID:?ERROR: run as a Slurm array}"


# input variables
tree_dir="$1"
vcf_dir="$2"
plink_bed_dir="$3"
pop_info_dir="$4"
anc_dir="$5"
global_anc_dir="$6"
sample_size="$7"
num_reps="$8"
msprime_model="$9"
chr="${10}"
genetic_map="${11}"
generation_time="${12}"
mutation_rate="${13}"
t_af_years="${14}"
t_ooa_years="${15}"
t_eu0_years="${16}"
t_eg_years="${17}"
r_eu0="${18}"
r_eu="${19}"
r_af="${20}"
n_a="${21}"
n_af1="${22}"
n_b="${23}"
n_eu0="${24}"
m_af_b="${25}"
m_af_eu="${26}"
admixture_time="${27}"
admix_generation_count="${28}"
admix_mixing_generation_count="${29}"
admix_ne_by_generation="${30}"
admix_afr_props_by_generation="${31}"
admix_eur_props_by_generation="${32}"
admix_prioradmix_props_by_generation="${33}"
admix_modern_growth_rate="${34}"
census_time_offset="${35}"
shift 35
shift 1 # skip the "--" from input arguments
pops=( "$@" )

# derived variables
rep="${SLURM_ARRAY_TASK_ID}"
prefix="${genetic_map}_${rep}_all"
tree_prefix="${tree_dir}/${prefix}"
vcf_path="${vcf_dir}/${prefix}.vcf"
plink_bed_prefix="${plink_bed_dir}/${prefix}"
pop_path="${pop_info_dir}/${genetic_map}_${rep}.pop"
if (( rep < 1 || rep > num_reps )); then
    echo "ERROR: invalid replicate ${rep}; expected 1..${num_reps}" >&2
    exit 1
fi

# run simulation using msprime via sim_3T.py
if [[ -s "${tree_prefix}.ts.tsz" && -s "${vcf_path}" \
    && -s "${plink_bed_prefix}.bed" ]]; then
    log_msg "simulation outputs exist for rep=${rep}; skipping"
    exit 0
fi

mkdir -p "${tree_dir}" "${vcf_dir}" "${plink_bed_dir}" "${pop_info_dir}" \
    "${anc_dir}" "${global_anc_dir}"
log_msg "running simulation/data generation for rep=${rep}"
python "${project_dir}/job_scripts/sim_model.py" \
    --tree-prefix "${tree_prefix}" \
    --vcf-path "${vcf_path}" \
    --pop-path "${pop_path}" \
    --anc-dir "${anc_dir}" \
    --global-anc-dir "${global_anc_dir}" \
    --sample-size "${sample_size}" \
    --seed "${rep}" \
    --msprime-model "${msprime_model}" \
    --chromosome "${chr}" \
    --genetic-map "${genetic_map}" \
    --generation-time "${generation_time}" \
    --mutation-rate "${mutation_rate}" \
    --t-af-years "${t_af_years}" \
    --t-ooa-years "${t_ooa_years}" \
    --t-eu0-years "${t_eu0_years}" \
    --t-eg-years "${t_eg_years}" \
    --r-eu0 "${r_eu0}" \
    --r-eu "${r_eu}" \
    --r-af "${r_af}" \
    --n-a "${n_a}" \
    --n-af1 "${n_af1}" \
    --n-b "${n_b}" \
    --n-eu0 "${n_eu0}" \
    --m-af-b "${m_af_b}" \
    --m-af-eu "${m_af_eu}" \
    --admixture-time "${admixture_time}" \
    --admix-generation-count "${admix_generation_count}" \
    --admix-mixing-generation-count "${admix_mixing_generation_count}" \
    --admix-ne-by-generation "${admix_ne_by_generation}" \
    --admix-afr-props-by-generation "${admix_afr_props_by_generation}" \
    --admix-eur-props-by-generation "${admix_eur_props_by_generation}" \
    --admix-prioradmix-props-by-generation \
    "${admix_prioradmix_props_by_generation}" \
    --admix-modern-growth-rate "${admix_modern_growth_rate}" \
    --census-time-offset "${census_time_offset}" \
    --pops "${pops[@]}"

# finalize files for downstream tasks
if [[ ! -s "${tree_prefix}.ts.tsz" ]]; then
    log_msg "compressing tree sequence for rep=${rep}"
    tszip "${tree_prefix}.ts"
fi
if [[ ! -s "${tree_prefix}.ts.tsz" ]]; then
    echo "ERROR: failed to create ${tree_prefix}.ts.tsz" >&2
    exit 1
fi
if [[ -s "${tree_prefix}.ts" ]]; then
    rm "${tree_prefix}.ts"
fi

log_msg "creating PLINK BED for rep=${rep}"
plink2 --vcf "${vcf_path}" --make-bed --out "${plink_bed_prefix}"

log_msg "done with simulation/data generation for rep=${rep}"

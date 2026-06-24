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
#SBATCH --output=/home1/karatas/logs/simOOANAA/simOOANAA.%A_%a.%x.out
#SBATCH --error=/home1/karatas/logs/simOOANAA/simOOANAA.%A_%a.%x.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=karatas@usc.edu


set -euo pipefail

module purge
ml gcc/13.3.0 htslib/1.19.1 bcftools/1.19 plink2/2.00a4.3 conda
source /apps/conda/miniforge3/25.3.0/etc/profile.d/conda.sh
conda activate OOA_NAAdmixture
export PATH="${HOME}/.conda/envs/OOA_NAAdmixture/bin:${PATH}"

project_dir="$(pwd)"
source "${project_dir}/other_scripts/log_msg.sh"

: "${SLURM_ARRAY_TASK_ID:?ERROR: run as a Slurm array}"


# input variables
tree_dir="$1"
pickled_demo_meta="$2"
vcf_dir="$3"
plink_bed_dir="$4"
pop_info_dir="$5"
anc_dir="$6"
global_anc_dir="$7"
sample_size="$8"
num_reps="$9"
msprime_model="${10}"
chr="${11}"
genetic_map="${12}"
generation_time="${13}"
mutation_rate="${14}"
t_af_years="${15}"
t_ooa_years="${16}"
t_eu0_years="${17}"
t_eg_years="${18}"
r_eu0="${19}"
r_eu="${20}"
r_af="${21}"
n_a="${22}"
n_af1="${23}"
n_b="${24}"
n_eu0="${25}"
m_af_b="${26}"
m_af_eu="${27}"
admixture_time="${28}"
admix_generation_count="${29}"
admix_mixing_generation_count="${30}"
admix_ne_by_generation="${31}"
admix_afr_props_by_generation="${32}"
admix_eur_props_by_generation="${33}"
admix_prioradmix_props_by_generation="${34}"
admix_modern_growth_rate="${35}"
census_time_offset="${36}"
shift 36
shift 1 # skip the "--" from input arguments
pops=( "$@" )

# derived variables
rep="${SLURM_ARRAY_TASK_ID}"
num_threads="${SLURM_CPUS_PER_TASK:-1}"
prefix="${genetic_map}_${rep}_all"
tree_prefix="${tree_dir}/${prefix}"
pickle_prefix="${pickled_demo_meta}/${prefix}"
vcf_path="${vcf_dir}/${prefix}.vcf"
vcf_gz_path="${vcf_path}.gz"
plink_vcf_path="${vcf_dir}/${prefix}.biallelic_snps.vcf.gz"
plink_bed_prefix="${plink_bed_dir}/${prefix}"
sample_metadata_path="${pop_info_dir}/${genetic_map}_${rep}.sample_metadata.tsv"
if (( rep < 1 || rep > num_reps )); then
    echo "ERROR: invalid replicate ${rep}; expected 1..${num_reps}" >&2
    exit 1
fi

# run simulation using msprime via sim_3T.py
if [[ -s "${tree_prefix}.ts.tsz" && -s "${vcf_gz_path}" \
    && -s "${plink_bed_prefix}.bed" && -s "${plink_bed_prefix}.bim" \
    && -s "${plink_bed_prefix}.fam" ]]; then
    log_msg "simulation outputs exist for rep=${rep}; skipping"
    exit 0
fi

mkdir -p "${tree_dir}" "${pickled_demo_meta}" "${vcf_dir}" \
    "${plink_bed_dir}" "${pop_info_dir}" "${anc_dir}" "${global_anc_dir}"
log_msg "running simulation/data generation for rep=${rep}"
python "${project_dir}/job_scripts/sim_model.py" \
    --tree-prefix "${tree_prefix}" \
    --pickle-prefix "${pickle_prefix}" \
    --vcf-path "${vcf_path}" \
    --sample-metadata-path "${sample_metadata_path}" \
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

log_msg "compressing and indexing all-sample VCF for rep=${rep}"
if [[ ! -s "${vcf_gz_path}" ]]; then
    bgzip -f "${vcf_path}"
fi
tabix -f -p vcf "${vcf_gz_path}"

log_msg "filtering biallelic SNPs for PLINK for rep=${rep}"
bcftools norm \
    --rm-dup all \
    --threads "${num_threads}" \
    "${vcf_gz_path}" |
    bcftools view \
        --types snps \
        --min-alleles 2 \
        --max-alleles 2 \
        --threads "${num_threads}" \
        -Oz \
        -o "${plink_vcf_path}"
tabix -p vcf "${plink_vcf_path}"

log_msg "creating PLINK BED for rep=${rep}"
plink2 \
    --vcf "${plink_vcf_path}" \
    --threads "${num_threads}" \
    --make-bed \
    --out "${plink_bed_prefix}"

if [[ ! -s "${plink_bed_prefix}.bed" || ! -s "${plink_bed_prefix}.bim" \
    || ! -s "${plink_bed_prefix}.fam" ]]; then
    echo "ERROR: failed to create complete PLINK BED set" >&2
    exit 1
fi

log_msg "done with simulation/data generation for rep=${rep}"

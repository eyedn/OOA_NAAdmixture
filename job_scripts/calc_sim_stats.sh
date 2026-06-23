#!/usr/bin/env bash

###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology 
#           Mooney Lab
#           ---
#           run_stats.sh
###############################################################################

#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --partition=qcb
#SBATCH --account=jazlynmo_738
#SBATCH --nodes=1
#SBATCH --output=/home1/karatas/logs/ooa_stats.%A_%a.%x.out
#SBATCH --error=/home1/karatas/logs/ooa_stats.%A_%a.%x.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=karatas@usc.edu


set -euo pipefail

module purge
ml gcc/13.3.0 htslib/1.19.1 conda
source /apps/conda/miniforge3/25.3.0/etc/profile.d/conda.sh
conda activate OOA_NAAdmixture
export PATH="${HOME}/.conda/envs/OOA_NAAdmixture/bin:${PATH}"

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
project_dir="$(cd "${script_dir}/.." && pwd)"
source "${project_dir}/other_scripts/log_msg.sh"

: "${SLURM_ARRAY_TASK_ID:?ERROR: run as a Slurm array}"


# input variables
tree_dir="$1"
plink_bed_dir="$2"
pop_info_dir="$3"
admixture_dir="$4"
king_dir="$5"
stats_dir="$6"
sample_size="$7"
num_reps="$8"
chr="$9"
genetic_map="${10}"
mutation_rate="${11}"
kin_cutoff="${12}"
shift 12
shift 1 # skip the "--" from input arguments
pops=( "$@" )

# derived variables
rep="${SLURM_ARRAY_TASK_ID}"
num_threads="${SLURM_CPUS_PER_TASK:-1}"
prefix="${genetic_map}_${rep}_all"
tree_tsz_path="${tree_dir}/${prefix}.ts.tsz"
plink_bed_prefix="${plink_bed_dir}/${prefix}"
pop_path="${pop_info_dir}/${genetic_map}_${rep}.pop"
admixture_prefix="${admixture_dir}/${prefix}"
if (( rep < 1 || rep > num_reps )); then
    echo "ERROR: invalid replicate ${rep}; expected 1..${num_reps}" >&2
    exit 1
fi

# run "admixture --supervised"
mkdir -p "${admixture_dir}" "${king_dir}" "${stats_dir}"
log_msg "running supervised ADMIXTURE for rep=${rep}"
cp "${pop_path}" "${admixture_prefix}.pop"
ln -sf "${plink_bed_prefix}.bed" "${admixture_prefix}.bed"
ln -sf "${plink_bed_prefix}.bim" "${admixture_prefix}.bim"
ln -sf "${plink_bed_prefix}.fam" "${admixture_prefix}.fam"
(
    cd "${admixture_dir}"
    "${HOME}/software/ADMIXTURE/admixture_linux-1.4.0/admixture" \
        --supervised \
        -j"${num_threads}" \
        -s "${rep}" \
        "${prefix}.bed" \
        2
)

# calculated KING kinships coefficients
log_msg "running within-population KING for rep=${rep}"
for pop in "${pops[@]}"; do
    subset_path="${king_dir}/${genetic_map}_${rep}_${pop}.subset"
    out_prefix="${king_dir}/${genetic_map}_${rep}_${pop}"
    python "${project_dir}/job_scripts/write_pop_subset.py" \
        --subset-path "${subset_path}" \
        --pop "${pop}" \
        --sample-size "${sample_size}" \
        --pops "${pops[@]}"
    plink2 \
        --bfile "${plink_bed_prefix}" \
        --keep "${subset_path}" \
        --threads "${num_threads}" \
        --king-cutoff "${kin_cutoff}" \
        --make-king \
        --make-king-table \
        --out "${out_prefix}"
done

# calculate summaries on pi, theta, sfs, ld, and kinship
log_msg "parsing stats for rep=${rep}"
python "${project_dir}/job_scripts/calc_sim_stats.py" \
    --rep "${rep}" \
    --tree-tsz-path "${tree_tsz_path}" \
    --plink-bed-prefix "${plink_bed_prefix}" \
    --admixture-dir "${admixture_dir}" \
    --king-dir "${king_dir}" \
    --stats-dir "${stats_dir}" \
    --sample-size "${sample_size}" \
    --chr "${chr}" \
    --genetic-map "${genetic_map}" \
    --mutation-rate "${mutation_rate}" \
    --pops "${pops[@]}"

log_msg "done with statistics for rep=${rep}"

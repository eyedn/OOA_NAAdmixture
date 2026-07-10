#!/usr/bin/env bash

###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology 
#           Mooney Lab
#           ---
#           combine_sim_stats.sh
###############################################################################


#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=1
#SBATCH --partition=qcb
#SBATCH --account=jazlynmo_738
#SBATCH --nodes=1
#SBATCH --output=/home1/karatas/logs/combOOANAA/combOOANAA.%A.%x.out
#SBATCH --error=/home1/karatas/logs/combOOANAA/combOOANAA.%A.%x.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=karatas@usc.edu


set -euo pipefail

module purge
ml gcc/13.3.0 htslib/1.19.1 conda
source /apps/conda/miniforge3/25.3.0/etc/profile.d/conda.sh
conda activate OOA_NAAdmixture
export PATH="${HOME}/.conda/envs/OOA_NAAdmixture/bin:${PATH}"

project_dir="$(pwd)"
source "${project_dir}/other_scripts/log_msg.sh"


# input variables
stats_dir="$1"
num_reps="$2"
shift 2
extra_args=( "$@" )

has_chromosomes=false
for extra_arg in "${extra_args[@]}"; do
    if [[ "${extra_arg}" == "--chromosomes" ]]; then
        has_chromosomes=true
        break
    fi
done

if [[ "${has_chromosomes}" == true && -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    extra_args+=( "--chrom-index" "${SLURM_ARRAY_TASK_ID}" )
fi


# combine per-replicate summaries
log_msg "combining simulation statistics across ${num_reps} replicates"
python "${project_dir}/job_scripts/python_utils/combine_sim_stats.py" \
    --stats-dir "${stats_dir}" \
    --num-reps "${num_reps}" \
    "${extra_args[@]}"

log_msg "done combining simulation statistics"

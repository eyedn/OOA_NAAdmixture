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

#SBATCH --time=4:00:00
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


# combine per-replicate summaries
log_msg "combining simulation statistics across ${num_reps} replicates"
python "${project_dir}/job_scripts/combine_sim_stats.py" \
    --stats-dir "${stats_dir}" \
    --num-reps "${num_reps}"

log_msg "done combining simulation statistics"

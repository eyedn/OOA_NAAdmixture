#!/usr/bin/env bash

###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           combine_sim_chr_stats.sh
###############################################################################

# workflow: combine simulation chromosome statistics across replicates.


##### set up ##################################################################
set -euo pipefail

# load required HPC modules and conda env
module purge
ml gcc/13.3.0 htslib/1.19.1 conda
source /apps/conda/miniforge3/25.3.0/etc/profile.d/conda.sh
conda activate OOA_NAAdmixture
export PATH="${HOME}/.conda/envs/OOA_NAAdmixture/bin:${PATH}"

# load shared functions; note, all scripts should exist in the execution repo
project_dir="$(pwd)"
source "${project_dir}/other_scripts/log_msg.sh"

# require execution as a Slurm array task
: "${SLURM_ARRAY_TASK_ID:?ERROR: run as a Slurm array}"


##### variables ###############################################################
# read the statistics directory, replicate count, and chromosome list.
stats_dir="$1"
num_reps="$2"
shift 2
shift 1 # skip the "--" from input arguments
chroms=( "$@" )


##### combine statistics ######################################################
# combine per-replicate summaries
log_msg "combining simulation chromosome statistics across replicates"
python "${project_dir}/job_scripts/python_utils/combine_sim_chr_stats.py" \
    --stats-dir "${stats_dir}" \
    --num-reps "${num_reps}" \
    --chrom-index "${SLURM_ARRAY_TASK_ID}" \
    --chroms "${chroms[@]}"

log_msg "done combining simulation chromosome statistics"

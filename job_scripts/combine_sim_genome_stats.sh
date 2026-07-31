#!/usr/bin/env bash

###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           combine_sim_genome_stats.sh
###############################################################################

# workflow: combine simulation genome statistics across replicates.


##### set up ##################################################################
set -euo pipefail

conda_env="$1"
shift

# load the software environment used by the table-combination helper.
module purge
ml gcc/13.3.0 htslib/1.19.1 conda
source /apps/conda/miniforge3/25.3.0/etc/profile.d/conda.sh
conda activate "${conda_env}"

# load shared timestamped logging.
project_dir="$(pwd)"
source "${project_dir}/other_scripts/log_msg.sh"


##### variables ###############################################################
# read the statistics directory and expected replicate count.
stats_dir="$1"
num_reps="$2"


##### combine statistics ######################################################
# combine per-replicate genome summaries into final TSV and Parquet tables.
log_msg "combining simulation genome statistics across replicates"
python "${project_dir}/job_scripts/python_utils/combine_sim_genome_stats.py" \
    --stats-dir "${stats_dir}" \
    --num-reps "${num_reps}"

log_msg "done combining simulation genome statistics"

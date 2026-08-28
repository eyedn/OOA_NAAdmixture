#!/usr/bin/env bash

###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           regenerate_onekg_genome_sfs.sh
###############################################################################

# workflow: aggregate projected chromosome SFS into canonical genome outputs.


##### set up ##################################################################
set -euo pipefail
conda_env="$1"
stats_dir="$2"
shift 2
shift
chroms=()
while [[ "$1" != "--" ]]; do
    chroms+=( "$1" )
    shift
done
shift
pops=( "$@" )

module purge
ml gcc/13.3.0 conda
source /apps/conda/miniforge3/25.3.0/etc/profile.d/conda.sh
conda activate "${conda_env}"
project_dir="$(pwd)"
source "${project_dir}/other_scripts/log_msg.sh"


##### aggregate SFS ###########################################################
log_msg "aggregating projected 1000 Genomes chromosome SFS"
python "${project_dir}/job_scripts/python_utils/regenerate_onekg_genome_sfs.py" \
    --stats-dir "${stats_dir}" \
    --chroms "${chroms[@]}" \
    --pops "${pops[@]}"
log_msg "completed projected 1000 Genomes genome SFS"

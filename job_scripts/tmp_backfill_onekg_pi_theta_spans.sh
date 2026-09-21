#!/usr/bin/env bash

###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           tmp_backfill_onekg_pi_theta_spans.sh
###############################################################################

# workflow: run the guarded post-hoc empirical genome pi/theta span repair.


##### set up ##################################################################
set -euo pipefail

conda_env="$1"
stats_dir="$2"
shift 2

module purge
ml gcc/13.3.0 plink2/2.00a4.3 conda
source /apps/conda/miniforge3/25.3.0/etc/profile.d/conda.sh
conda activate "${conda_env}"


##### arguments ###############################################################
if [[ "$1" != "--" ]]; then
    echo "ERROR: expected chromosome delimiter" >&2
    exit 1
fi
shift
chroms=()
while [[ "$1" != "--" ]]; do
    chroms+=( "$1" )
    shift
done
shift
pops=( "$@" )
if (( ${#chroms[@]} == 0 || ${#pops[@]} == 0 )); then
    echo "ERROR: at least one chromosome and population are required" >&2
    exit 1
fi


##### backfill ###############################################################
project_dir="$(pwd)"
python "${project_dir}/job_scripts/python_utils/"\
"tmp_backfill_onekg_pi_theta_spans.py" \
    --stats-dir "${stats_dir}" \
    --chroms "${chroms[@]}" \
    --pops "${pops[@]}"

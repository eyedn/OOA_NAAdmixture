#!/usr/bin/env bash

###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           regenerate_sim_sfs.sh
###############################################################################

# workflow: regenerate one projected simulation chromosome-replicate SFS.


##### set up ##################################################################
set -euo pipefail

conda_env="$1"
tree_dir="$2"
stats_dir="$3"
sample_size="$4"
sfs_size="$5"
num_reps="$6"
genetic_map="$7"
shift 7
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
source "${project_dir}/other_scripts/map_slurm_task.sh"
: "${SLURM_ARRAY_TASK_ID:?ERROR: run as a Slurm array}"


##### regenerate SFS ##########################################################
read -r chr rep < <(
    map_slurm_task_to_chr_rep \
        "${SLURM_ARRAY_TASK_ID}" "${num_reps}" -- "${chroms[@]}"
)
prefix="${genetic_map}_${rep}_chr${chr}_all"
log_msg "regenerating simulation SFS rep=${rep} chr=${chr}"
python "${project_dir}/job_scripts/python_utils/regenerate_sim_sfs.py" \
    --rep "${rep}" \
    --chrom "${chr}" \
    --tree-tsz-path "${tree_dir}/${prefix}.ts.tsz" \
    --stats-dir "${stats_dir}" \
    --sample-size "${sample_size}" \
    --sfs-size "${sfs_size}" \
    --pops "${pops[@]}"
log_msg "completed simulation SFS rep=${rep} chr=${chr}"

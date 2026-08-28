#!/usr/bin/env bash

###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           tmp_submit_regenerate_sim_sfs.sh
###############################################################################

# workflow: regenerate only canonical projected simulation SFS outputs.


##### set up ##################################################################
set -euo pipefail
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${script_dir}/other_scripts/const.sh"
source "${script_dir}/other_scripts/log_msg.sh"
mkdir -p "${STATS_DIR}"


##### chromosome-replicate recovery ##########################################
array_size=$((NUM_REPS * ${#CHROMS[@]}))
stats_jname="regenSimSFS"
mkdir -p "/home1/karatas/logs/${stats_jname}"
stats_jid=$(sbatch \
    --parsable \
    --chdir="${script_dir}" \
    --job-name="${stats_jname}" \
    --array="1-${array_size}%${MAX_JOBS}" \
    --cpus-per-task="${SIM_STATS_CPUS_PER_TASK}" \
    --mem="${SIM_STATS_MEM}" \
    --time=1-00:00:00 \
    --partition="${PARTITION}" \
    --account="${ACCOUNT}" \
    --nodes=1 \
    --output="/home1/karatas/logs/${stats_jname}/%A_%a.%x.out" \
    --error="/home1/karatas/logs/${stats_jname}/%A_%a.%x.err" \
    --mail-type="${MAIL_TYPE}" \
    --mail-user="${MAIL_USER}" \
    "job_scripts/regenerate_sim_sfs.sh" \
        "${OOA_NAADMIXTURE_CONDA}" \
        "${TREE_DIR}" \
        "${STATS_DIR}" \
        "${SAMPLE_SIZE}" \
        "${SFS_SIZE}" \
        "${NUM_REPS}" \
        "${GENETIC_MAP}" \
        -- "${CHROMS[@]}" -- "${POPS[@]}"
)


##### canonical chromosome combination #######################################
comb_jname="combRegenSimSFS"
mkdir -p "/home1/karatas/logs/${comb_jname}"
comb_jid=$(sbatch \
    --parsable \
    --dependency="afterok:${stats_jid}" \
    --chdir="${script_dir}" \
    --job-name="${comb_jname}" \
    --array="1-${#CHROMS[@]}%${MAX_JOBS}" \
    --cpus-per-task=1 \
    --mem="${SIM_COMB_MEM}" \
    --time=1-00:00:00 \
    --partition="${PARTITION}" \
    --account="${ACCOUNT}" \
    --nodes=1 \
    --output="/home1/karatas/logs/${comb_jname}/%A_%a.%x.out" \
    --error="/home1/karatas/logs/${comb_jname}/%A_%a.%x.err" \
    --mail-type="${MAIL_TYPE}" \
    --mail-user="${MAIL_USER}" \
    "job_scripts/combine_sim_chr_stats.sh" \
        "${OOA_NAADMIXTURE_CONDA}" \
        "${STATS_DIR}" \
        "${NUM_REPS}" \
        -- "${CHROMS[@]}" -- "sfs"
)
log_msg "submitted simulation SFS recovery; stats=${stats_jid} combine=${comb_jid}"

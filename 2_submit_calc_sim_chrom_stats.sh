#!/usr/bin/env bash

###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           2_submit_calc_sim_chrom_stats.sh
###############################################################################

# workflow: submit simulation chromosome-statistics jobs after simulation.


##### set up ##################################################################
set -euo pipefail

# determine repo location; note, all scripts should exist in the execution repo
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# source shared constants and functions
source "${script_dir}/other_scripts/const.sh"
source "${script_dir}/other_scripts/log_msg.sh"


##### chromosome statistics jobs ##############################################
# create output directories for ADMIXTURE, KING, and stats outputs
mkdir -p "${ADMIXTURE_DIR}" "${KING_DIR}" "${STATS_DIR}"

# pass simulation inputs and analysis parameters to worker script; array loops
# over chromosome-by-replicate combinations
log_msg "submitting OOA_NAAdmixture chromosome statistics array
OUT_DIR=${OUT_DIR}
num_reps=${NUM_REPS}
sample_size=${SAMPLE_SIZE}
chroms=${CHROMS[*]}
pops=${POPS[*]}"
array_size=$((NUM_REPS * ${#CHROMS[@]}))
stats_jname="calcChrOOANAA"
mkdir -p "/home1/karatas/logs/${stats_jname}"
stats_jid=$(sbatch \
    --parsable \
    --chdir="$script_dir" \
    --job-name="$stats_jname" \
    --array="1-${array_size}%${MAX_JOBS}" \
    --cpus-per-task="${STATS_CPUS_PER_TASK}" \
    --mem="${STATS_MEM}" \
    --time=1-00:00:00 \
    --partition="${PARTITION}" \
    --account="${ACCOUNT}" \
    --nodes=1 \
    --output="/home1/karatas/logs/${stats_jname}/%A_%a.%x.out" \
    --error="/home1/karatas/logs/${stats_jname}/%A_%a.%x.err" \
    --mail-type="${MAIL_TYPE}" \
    --mail-user="${MAIL_USER}" \
    "job_scripts/calc_sim_chr_rep_stats.sh" \
        "${TREE_DIR}" \
        "${PLINK_BED_DIR}" \
        "${POP_INFO_DIR}" \
        "${ADMIXTURE_DIR}" \
        "${GLOBAL_ANC_DIR}" \
        "${KING_DIR}" \
        "${STATS_DIR}" \
        "${SAMPLE_SIZE}" \
        "${NUM_REPS}" \
        "${GENETIC_MAP}" \
        "${MUTATION_RATE}" \
        "${KIN_CUTOFF=}" \
        "${ADMIXTURE_LD_WINDOW}" \
        "${ADMIXTURE_LD_STEP}" \
        "${ADMIXTURE_LD_R2}" \
        -- \
        "${CHROMS[@]}" \
        -- \
        "${POPS[@]}"
)
log_msg "submitted chromosome statistics array; jid=${stats_jid}"


##### chromosome statistics combination #######################################
# pass stats outputs worker script that combines the results of the previous
# job; array loops through chrom
comb_jname="combChrOOANAA"
mkdir -p "/home1/karatas/logs/${comb_jname}"
comb_jid=$(sbatch \
    --parsable \
    --dependency="afterok:${stats_jid}" \
    --chdir="${script_dir}" \
    --job-name="$comb_jname" \
    --array="1-${#CHROMS[@]}%${MAX_JOBS}" \
    --mem="${COMB_MEM}" \
    --time=1-00:00:00 \
    --cpus-per-task=1 \
    --partition="${PARTITION}" \
    --account="${ACCOUNT}" \
    --nodes=1 \
    --output="/home1/karatas/logs/${comb_jname}/%A_%a.%x.out" \
    --error="/home1/karatas/logs/${comb_jname}/%A_%a.%x.err" \
    --mail-type="${MAIL_TYPE}" \
    --mail-user="${MAIL_USER}" \
    "job_scripts/combine_sim_chr_stats.sh" \
        "${STATS_DIR}" \
        "${NUM_REPS}" \
        -- \
        "${CHROMS[@]}"
)
log_msg "submitted chromosome statistics combine job; jid=${comb_jid}"

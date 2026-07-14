#!/usr/bin/env bash

###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           temp_submit_calc_genome_kinship_unrelated.sh
###############################################################################


set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${script_dir}/other_scripts/const.sh"
source "${script_dir}/other_scripts/log_msg.sh"


# submit the focused genome backfill array
mkdir -p "${KING_DIR}" "${STATS_DIR}"
stats_jname="tempGenomeKingUnrel"
mkdir -p "/home1/karatas/logs/${stats_jname}"
stats_jid=$(sbatch \
    --parsable \
    --chdir="${script_dir}" \
    --job-name="${stats_jname}" \
    --array="1-${NUM_REPS}%${MAX_JOBS}" \
    --cpus-per-task="${STATS_CPUS_PER_TASK}" \
    --mem="${STATS_MEM}" \
    --time=1-00:00:00 \
    --partition=qcb \
    --account=jazlynmo_738 \
    --nodes=1 \
    --output="/home1/karatas/logs/${stats_jname}/%A_%a.%x.out" \
    --error="/home1/karatas/logs/${stats_jname}/%A_%a.%x.err" \
    --mail-type=ALL \
    --mail-user=karatas@usc.edu \
    "job_scripts/temp_calc_genome_kinship_unrelated.sh" \
        "${PLINK_BED_DIR}" \
        "${KING_DIR}" \
        "${STATS_DIR}" \
        "${NUM_REPS}" \
        "${GENETIC_MAP}" \
        -- \
        "${POPS[@]}"
)

log_msg "submitted genome unrelated KING array; jid=${stats_jid}"

# combine only new genome tables after the full array succeeds
comb_jname="tempCombGenomeKingUnrel"
mkdir -p "/home1/karatas/logs/${comb_jname}"
comb_jid=$(sbatch \
    --parsable \
    --dependency="afterok:${stats_jid}" \
    --chdir="${script_dir}" \
    --job-name="${comb_jname}" \
    --mem="${COMB_MEM}" \
    --time=1-00:00:00 \
    --cpus-per-task=1 \
    --partition=qcb \
    --account=jazlynmo_738 \
    --nodes=1 \
    --output="/home1/karatas/logs/${comb_jname}/%A.%x.out" \
    --error="/home1/karatas/logs/${comb_jname}/%A.%x.err" \
    --mail-type=ALL \
    --mail-user=karatas@usc.edu \
    "job_scripts/combine_sim_stats.sh" \
        "${STATS_DIR}" \
        "${NUM_REPS}" \
        --tables kinship_unrelated
)

log_msg "submitted genome unrelated KING combine; jid=${comb_jid}"

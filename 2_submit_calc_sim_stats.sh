#!/usr/bin/env bash

###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology 
#           Mooney Lab
#           ---
#           2_submit_calc_sim_stats.sh
###############################################################################


set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${script_dir}/other_scripts/const.sh"
source "${script_dir}/other_scripts/log_msg.sh"


# generate summaries of pi, watterson's theta, 1D/2D sfs, ld decay, and kinship
mkdir -p "${ADMIXTURE_DIR}" "${KING_DIR}" "${STATS_DIR}"
log_msg "submitting OOA_NAAdmixture statistics array
outdir=${OUTDIR}
num_reps=${NUM_REPS}
sample_size=${SAMPLE_SIZE}
pops=${POPS[*]}"

stats_jid=$(sbatch \
    --parsable \
    --chdir="${script_dir}" \
    --job-name="ooa_stats" \
    --array="1-${NUM_REPS}%${MAX_JOBS}" \
    "job_scripts/calc_sim_stats.sh" \
        "${TREE_DIR}" \
        "${PLINK_BED_DIR}" \
        "${POP_INFO_DIR}" \
        "${ADMIXTURE_DIR}" \
        "${KING_DIR}" \
        "${STATS_DIR}" \
        "${SAMPLE_SIZE}" \
        "${NUM_REPS}" \
        "${CHR}" \
        "${GENETIC_MAP}" \
        "${MUTATION_RATE}" \
        "${KIN_CUTOFF=}" \
        -- \
        "${POPS[@]}"
)

log_msg "submitted simulation statistics array; jid=${stats_jid}"

combine_jid=$(sbatch \
    --parsable \
    --dependency="afterok:${stats_jid}" \
    --chdir="${script_dir}" \
    --job-name="ooa_combine_stats" \
    "job_scripts/combine_sim_stats.sh" \
        "${STATS_DIR}" \
        "${NUM_REPS}"
)

log_msg "submitted simulation statistics combine job; jid=${combine_jid}"

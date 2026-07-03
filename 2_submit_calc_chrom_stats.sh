#!/usr/bin/env bash

###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology 
#           Mooney Lab
#           ---
#           2_submit_calc_chrom_stats.sh
###############################################################################


set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${script_dir}/other_scripts/const.sh"
source "${script_dir}/other_scripts/log_msg.sh"


# generate chromosome-level pi, theta, SFS, LD, KING, and ADMIXTURE summaries
mkdir -p "${ADMIXTURE_DIR}" "${KING_DIR}" "${STATS_DIR}"
array_size=$((NUM_REPS * ${#CHROMS[@]}))
log_msg "submitting OOA_NAAdmixture chromosome statistics array
outdir=${OUTDIR}
num_reps=${NUM_REPS}
sample_size=${SAMPLE_SIZE}
chroms=${CHROMS[*]}
pops=${POPS[*]}"

stats_jid=$(sbatch \
    --parsable \
    --chdir="${script_dir}" \
    --job-name="calcChrOOANAA" \
    --array="1-${array_size}%${MAX_JOBS}" \
    --cpus-per-task="${STATS_CPUS_PER_TASK}" \
    --mem="${STATS_MEM}" \
    "job_scripts/calc_chrom_stats.sh" \
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

chrom_combine_jid=$(sbatch \
    --parsable \
    --dependency="afterok:${stats_jid}" \
    --chdir="${script_dir}" \
    --job-name="combChrOOANAA" \
    --array="1-${#CHROMS[@]}%${MAX_JOBS}" \
    --mem="${COMB_MEM}" \
    "job_scripts/combine_sim_stats.sh" \
        "${STATS_DIR}" \
        "${NUM_REPS}" \
        --chromosomes \
        --chroms \
        "${CHROMS[@]}"
)

log_msg "submitted chromosome statistics combine job; jid=${chrom_combine_jid}"

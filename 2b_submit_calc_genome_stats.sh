#!/usr/bin/env bash

###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology 
#           Mooney Lab
#           ---
#           2b_submit_calc_genome_stats.sh
###############################################################################


set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${script_dir}/other_scripts/const.sh"
source "${script_dir}/other_scripts/log_msg.sh"


# aggregate chromosome statistics and run genome-level ADMIXTURE/KING
mkdir -p "${ADMIXTURE_DIR}" "${KING_DIR}" "${STATS_DIR}"
log_msg "submitting OOA_NAAdmixture genome statistics array
outdir=${OUTDIR}
num_reps=${NUM_REPS}
sample_size=${SAMPLE_SIZE}
chroms=${CHROMS[*]}
pops=${POPS[*]}"

genome_jid=$(sbatch \
    --parsable \
    --chdir="${script_dir}" \
    --job-name="calcGenomeOOANAA" \
    --array="1-${NUM_REPS}%${MAX_JOBS}" \
    --cpus-per-task="${STATS_CPUS_PER_TASK}" \
    --mem="${STATS_MEM}" \
    "job_scripts/calc_genome_stats.sh" \
        "${VCF_DIR}" \
        "${PLINK_BED_DIR}" \
        "${POP_INFO_DIR}" \
        "${ADMIXTURE_DIR}" \
        "${KING_DIR}" \
        "${STATS_DIR}" \
        "${SAMPLE_SIZE}" \
        "${NUM_REPS}" \
        "${GENETIC_MAP}" \
        "${KIN_CUTOFF=}" \
        "${ADMIXTURE_LD_WINDOW}" \
        "${ADMIXTURE_LD_STEP}" \
        "${ADMIXTURE_LD_R2}" \
        "${CHROMS[@]}" \
        -- \
        "${POPS[@]}"
)

log_msg "submitted genome statistics array; jid=${genome_jid}"

combine_jid=$(sbatch \
    --parsable \
    --dependency="afterok:${genome_jid}" \
    --chdir="${script_dir}" \
    --job-name="combOOANAA" \
    --mem="${COMB_MEM}" \
    "job_scripts/combine_sim_stats.sh" \
        "${STATS_DIR}" \
        "${NUM_REPS}"
)

log_msg "submitted genome statistics combine job; jid=${combine_jid}"

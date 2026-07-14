#!/usr/bin/env bash

###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology 
#           Mooney Lab
#           ---
#           3_submit_calc_genome_stats.sh
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

stats_jname="calcGenomeOOANAA"
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
        -- \
        "${CHROMS[@]}" \
        -- \
        "${POPS[@]}"
)

log_msg "submitted genome statistics array; jid=${stats_jid}"

comb_jname="combGenomeOOANAA"
mkdir -p "/home1/karatas/logs/${comb_jname}"
comb_jid=$(sbatch \
    --parsable \
    --chdir="${script_dir}" \
    --dependency="afterok:${stats_jid}" \
    --job-name="combGenomeOOANAA" \
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
        "${NUM_REPS}"
)

log_msg "submitted genome statistics combine job; jid=${comb_jid}"

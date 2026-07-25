#!/usr/bin/env bash

###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           3_submit_calc_sim_genome_stats.sh
###############################################################################

# workflow: submit per-replicate genome statistics after chromosome jobs, then
# submit the dependent combine job that writes genome-level summary tables.


##### set up ##################################################################
set -euo pipefail

# determine repo location; all referenced scripts run from this checkout.
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# source shared pipeline constants and timestamped logging.
source "${script_dir}/other_scripts/const.sh"
source "${script_dir}/other_scripts/log_msg.sh"

if [[ -z ${OOA_NAADMIXTURE_CONDA} ]]; then
    echo "ERROR: OOA_NAADMIXTURE_CONDA must be configured" >&2
    exit 1
fi
if [[ ! -x "${ADMIXTURE_EXEC}" ]]; then
    echo "ERROR: ADMIXTURE executable is not executable: ${ADMIXTURE_EXEC}" >&2
    exit 1
fi


##### genome statistics jobs ##################################################
# create shared output directories before submitting replicate-level workers.
mkdir -p "${ADMIXTURE_DIR}" "${KING_DIR}" "${STATS_DIR}"

# submit one worker per replicate; each worker merges all chromosome inputs.
log_msg "submitting OOA_NAAdmixture genome statistics array
OUT_DIR=${OUT_DIR}
num_reps=${NUM_REPS}
sample_size=${SAMPLE_SIZE}
chroms=${CHROMS[*]}
pops=${POPS[*]}"

stats_jname="calcGenomeOOANAA"
mkdir -p "/home1/karatas/logs/${stats_jname}"

# pass chromosome artifacts and analysis settings to the genome worker.
stats_jid=$(sbatch \
    --parsable \
    --chdir="${script_dir}" \
    --job-name="${stats_jname}" \
    --array="1-${NUM_REPS}%${MAX_JOBS}" \
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
    "job_scripts/calc_sim_genome_rep_stats.sh" \
        "${OOA_NAADMIXTURE_CONDA}" \
        "${ADMIXTURE_EXEC}" \
        "${VCF_DIR}" \
        "${PLINK_BED_DIR}" \
        "${POP_INFO_DIR}" \
        "${ADMIXTURE_DIR}" \
        "${KING_DIR}" \
        "${STATS_DIR}" \
        "${SAMPLE_SIZE}" \
        "${NUM_REPS}" \
        "${GENETIC_MAP}" \
        "${KIN_CUTOFF}" \
        "${ADMIXTURE_LD_WINDOW}" \
        "${ADMIXTURE_LD_STEP}" \
        "${ADMIXTURE_LD_R2}" \
        -- \
        "${CHROMS[@]}" \
        -- \
        "${POPS[@]}"
)
log_msg "submitted genome statistics array; jid=${stats_jid}"


##### genome statistics combination ###########################################
# combine per-replicate genome tables only after all workers succeed.
comb_jname="combGenomeOOANAA"
mkdir -p "/home1/karatas/logs/${comb_jname}"
comb_jid=$(sbatch \
    --parsable \
    --chdir="${script_dir}" \
    --dependency="afterok:${stats_jid}" \
    --job-name="combGenomeOOANAA" \
    --mem="${SIM_COMB_MEM}" \
    --time=1-00:00:00 \
    --cpus-per-task=1 \
    --partition="${PARTITION}" \
    --account="${ACCOUNT}" \
    --nodes=1 \
    --output="/home1/karatas/logs/${comb_jname}/%A.%x.out" \
    --error="/home1/karatas/logs/${comb_jname}/%A.%x.err" \
    --mail-type="${MAIL_TYPE}" \
    --mail-user="${MAIL_USER}" \
    "job_scripts/combine_sim_genome_stats.sh" \
        "${OOA_NAADMIXTURE_CONDA}" \
        "${STATS_DIR}" \
        "${NUM_REPS}"
)
log_msg "submitted genome statistics combine job; jid=${comb_jid}"

#!/usr/bin/env bash

###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           5_submit_calc_onekg_genome_stats.sh
###############################################################################

# workflow: validate chromosome outputs, calculate per-population genome
# statistics, then combine the resulting empirical genome tables.


##### set up ##################################################################
set -euo pipefail

# determine repo location; all referenced scripts run from this checkout.
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# source shared empirical constants and timestamped logging.
source "${script_dir}/other_scripts/const.sh"
source "${script_dir}/other_scripts/log_msg.sh"


##### genome statistics jobs ##################################################
# submit one population-specific worker after all chromosome handoffs validate.
pop_jname="statsOnekgGenomePop"
mkdir -p "/home1/karatas/logs/${pop_jname}"
pop_jid=$(sbatch \
    --parsable \
    --chdir="${script_dir}" \
    --job-name="${pop_jname}" \
    --array="1-${#ONEKG_POPS[@]}%${MAX_JOBS}" \
    --cpus-per-task="${STATS_CPUS_PER_TASK}" \
    --mem="${STATS_MEM}" \
    --time=1-00:00:00 \
    --partition="${PARTITION}" \
    --account="${ACCOUNT}" \
    --nodes=1 \
    --output="/home1/karatas/logs/${pop_jname}/%A_%a.%x.out" \
    --error="/home1/karatas/logs/${pop_jname}/%A_%a.%x.err" \
    --mail-type="${MAIL_TYPE}" \
    --mail-user="${MAIL_USER}" \
    "job_scripts/calc_onekg_genome_pop_stats.sh" \
        "${ONEKG_OUT_PLINK_BED_DIR}" \
        "${ONEKG_OUT_POP_INFO_DIR}" \
        "${ONEKG_OUT_KING_DIR}" \
        "${ONEKG_OUT_STATS_DIR}" \
        -- \
        "${CHROMS[@]}" \
        -- \
        "${ONEKG_POPS[@]}"
)
log_msg "submitted 1000 Genomes genome stats. job; jid=${pop_jid}"


##### genome statistics combination ###########################################
# combine all population results only after every genome worker succeeds.
comb_jname="combOnekgGenome"
mkdir -p "/home1/karatas/logs/${comb_jname}"
comb_jid=$(sbatch \
    --parsable \
    --dependency="afterok:${pop_jid}" \
    --chdir="${script_dir}" \
    --job-name="${comb_jname}" \
    --mem="${COMB_MEM}" \
    --time=1-00:00:00 \
    --partition="${PARTITION}" \
    --account="${ACCOUNT}" \
    --nodes=1 \
    --output="/home1/karatas/logs/${comb_jname}/%A.%x.out" \
    --error="/home1/karatas/logs/${comb_jname}/%A.%x.err" \
    --mail-type="${MAIL_TYPE}" \
    --mail-user="${MAIL_USER}" \
    "job_scripts/combine_onekg_genome_stats.sh" \
        "${ONEKG_UNRELS_FILE}" \
        "${ONEKG_PLINK_FAM_FILE}" \
        "${ONEKG_OUT_PLINK_BED_DIR}" \
        "${ONEKG_OUT_ADMIXTURE_DIR}" \
        "${ONEKG_OUT_STATS_DIR}" \
        "${ADMIXTURE_LD_WINDOW}" \
        "${ADMIXTURE_LD_STEP}" \
        "${ADMIXTURE_LD_R2}" \
        -- \
        "${CHROMS[@]}" \
        -- \
        "${ONEKG_POPS[@]}"
)
log_msg "submitted 1000 Genomes genome comb. job; jid=${comb_jid}"

#!/usr/bin/env bash

###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           tmp_submit_regenerate_onekg_sfs.sh
###############################################################################

# workflow: regenerate only canonical empirical chromosome and genome SFS.


##### set up ##################################################################
set -euo pipefail
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${script_dir}/other_scripts/const.sh"
source "${script_dir}/other_scripts/log_msg.sh"
mkdir -p "${ONEKG_OUT_STATS_DIR}"


##### chromosome recovery #####################################################
stats_jname="regen1kGSFS"
mkdir -p "/home1/karatas/logs/${stats_jname}"
stats_jid=$(sbatch \
    --parsable \
    --chdir="${script_dir}" \
    --job-name="${stats_jname}" \
    --array="1-${#CHROMS[@]}%${MAX_JOBS}" \
    --cpus-per-task="${ONEKG_STATS_CPUS_PER_TASK}" \
    --mem="${ONEKG_STATS_MEM}" \
    --time=1-00:00:00 \
    --partition="${PARTITION}" \
    --account="${ACCOUNT}" \
    --nodes=1 \
    --output="/home1/karatas/logs/${stats_jname}/%A_%a.%x.out" \
    --error="/home1/karatas/logs/${stats_jname}/%A_%a.%x.err" \
    --mail-type="${MAIL_TYPE}" \
    --mail-user="${MAIL_USER}" \
    "job_scripts/regenerate_onekg_chr_sfs.sh" \
        "${OOA_NAADMIXTURE_CONDA}" \
        "${ONEKG_VCF_PREFIX}" \
        "${ONEKG_VCF_SUFFIX}" \
        "${ONEKG_UNRELS_FILE}" \
        "${ONEKG_PLINK_FAM_FILE}" \
        "${ONEKG_OUT_STATS_DIR}" \
        "${SFS_SIZE}" \
        "${SFS_SIZE_POP_REF}" \
        -- "${CHROMS[@]}" -- "${ONEKG_POPS[@]}"
)


##### dependent genome aggregation ###########################################
genome_jname="regen1kGGenomeSFS"
mkdir -p "/home1/karatas/logs/${genome_jname}"
genome_jid=$(sbatch \
    --parsable \
    --dependency="afterok:${stats_jid}" \
    --chdir="${script_dir}" \
    --job-name="${genome_jname}" \
    --cpus-per-task=1 \
    --mem="${ONEKG_COMB_MEM}" \
    --time=1-00:00:00 \
    --partition="${PARTITION}" \
    --account="${ACCOUNT}" \
    --nodes=1 \
    --output="/home1/karatas/logs/${genome_jname}/%A.%x.out" \
    --error="/home1/karatas/logs/${genome_jname}/%A.%x.err" \
    --mail-type="${MAIL_TYPE}" \
    --mail-user="${MAIL_USER}" \
    "job_scripts/regenerate_onekg_genome_sfs.sh" \
        "${OOA_NAADMIXTURE_CONDA}" \
        "${ONEKG_OUT_STATS_DIR}" \
        -- "${CHROMS[@]}" -- "${ONEKG_POPS[@]}"
)
log_msg "submitted empirical SFS recovery; stats=${stats_jid} genome=${genome_jid}"

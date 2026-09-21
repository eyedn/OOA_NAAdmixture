#!/usr/bin/env bash

###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           tmp_submit_backfill_onekg_pi_theta_spans.sh
###############################################################################

# workflow: submit one guarded post-hoc empirical genome pi/theta span repair.


##### set up ##################################################################
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${script_dir}/other_scripts/const.sh"


##### post-hoc backfill job ##################################################
job_name="backfillOnekgPiThetaSpan"
log_dir="/home1/karatas/logs/${job_name}"
mkdir -p "${log_dir}"
job_id=$(sbatch \
    --parsable \
    --chdir="${script_dir}" \
    --job-name="${job_name}" \
    --cpus-per-task="1" \
    --mem="${ONEKG_COMB_MEM}" \
    --time="1-00:00:00" \
    --partition="${PARTITION}" \
    --account="${ACCOUNT}" \
    --nodes="1" \
    --output="${log_dir}/%j.%x.out" \
    --error="${log_dir}/%j.%x.err" \
    --mail-type="${MAIL_TYPE}" \
    --mail-user="${MAIL_USER}" \
    "${script_dir}/job_scripts/tmp_backfill_onekg_pi_theta_spans.sh" \
        "${OOA_NAADMIXTURE_CONDA}" \
        "${ONEKG_OUT_STATS_DIR}" \
        -- \
        "${CHROMS[@]}" \
        -- \
        "${ONEKG_POPS[@]}"
)
printf 'submitted empirical genome pi/theta span backfill; jid=%s\n' "${job_id}"

#!/usr/bin/env bash

###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           sim_select_stdpopsim_models.sh
###############################################################################

# workflow: submit one chrom. 22 simulation for each selected stdpopsim model.


##### set up ##################################################################
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
project_dir="$(cd "${script_dir}/.." && pwd)"
source "${script_dir}/const.sh"
source "${script_dir}/log_msg.sh"

: "${OOA_NAADMIXTURE_CONDA:?ERROR: conda environment is not configured}"

select_out_dir="/home1/karatas/scratch/OOA_NAAdmixture_selectStspopsimModels"
log_dir="/home1/karatas/logs/selectStdpopsimModels"
mkdir -p "${log_dir}"


##### submit array ############################################################
log_msg "submitting selected stdpopsim models; output=${select_out_dir}"
sim_jid=$(sbatch \
    --parsable \
    --chdir="${project_dir}" \
    --job-name="selectStdpopsimModels" \
    --array="1-3%${MAX_JOBS}" \
    --cpus-per-task="${SIM_CPUS_PER_TASK}" \
    --mem="${SIM_MEM}" \
    --time=1-00:00:00 \
    --partition="${PARTITION}" \
    --account="${ACCOUNT}" \
    --nodes=1 \
    --output="${log_dir}/%A_%a.%x.out" \
    --error="${log_dir}/%A_%a.%x.err" \
    --mail-type="${MAIL_TYPE}" \
    --mail-user="${MAIL_USER}" \
    "${script_dir}/sim_select_stdpopsim_models_worker.sh" \
    "${OOA_NAADMIXTURE_CONDA}" \
    "${select_out_dir}" \
    "${GENETIC_MAP}" \
    "${MSPRIME_MODEL}"
)
log_msg "submitted selected stdpopsim model array; jid=${sim_jid}"

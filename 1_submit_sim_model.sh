#!/usr/bin/env bash

###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           1_submit_sim_model.sh
###############################################################################

# workflow: submit the chromosome-by-replicate simulation array and hand off
# generated tree sequences, VCFs, PLINK files, and ancestry metadata.


##### set up ##################################################################
set -euo pipefail

# determine repo location; note, all scripts should exist in the execution repo
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# source shared constants and functions
source "${script_dir}/other_scripts/const.sh"
source "${script_dir}/other_scripts/log_msg.sh"


##### pre-simulation prep #####################################################
# create comma-separated lists of admixture variables from constants for
# passing to worker scripts
join_by_comma() {
    local IFS=","
    printf "%s" "$*"
}
admix_ne_by_generation=$(join_by_comma "${ADMIX_NE_BY_GENERATION[@]}")
admix_afr_props_by_generation=$(
    join_by_comma "${ADMIX_AFR_PROPS_BY_GENERATION[@]}"
)
admix_eur_props_by_generation=$(
    join_by_comma "${ADMIX_EUR_PROPS_BY_GENERATION[@]}"
)
admix_prioradmix_props_by_generation=$(
    join_by_comma "${ADMIX_PRIORADMIX_PROPS_BY_GENERATION[@]}"
)


##### simulation ##############################################################
# create output directories for simulated tree sequences, msprime metadata,
# vcf, PLINK binaries, population metadata, local ancestry, global ancestry
mkdir -p "${TREE_DIR}" "${PICKLED_DEMO_META}" "${VCF_DIR}" \
    "${PLINK_BED_DIR}" "${POP_INFO_DIR}" "${ANC_DIR}" "${GLOBAL_ANC_DIR}"

# pass simulation and demographic parameters to simulation worker script
log_msg "submitting OOA_NAAdmixture simulation array
OUT_DIR=${OUT_DIR}
num_reps=${NUM_REPS}
sample_size=${SAMPLE_SIZE}
chroms=${CHROMS[*]}
genetic_map=${GENETIC_MAP}
pops=${POPS[*]}"
array_size=$((NUM_REPS * ${#CHROMS[@]}))
jname="simOOANAA"
mkdir -p "/home1/karatas/logs/${jname}"
sim_jid=$(sbatch \
    --parsable \
    --chdir="${script_dir}" \
    --job-name="${jname}" \
    --array="1-${array_size}%${MAX_JOBS}" \
    --cpus-per-task="${SIM_CPUS_PER_TASK}" \
    --mem="${SIM_MEM}" \
    --time=1-00:00:00 \
    --partition=qcb \
    --account=jazlynmo_738 \
    --nodes=1 \
    --output="/home1/karatas/logs/${jname}/%A_%a.%x.out" \
    --error="/home1/karatas/logs/${jname}/%A_%a.%x.err" \
    --mail-type=ALL \
    --mail-user=karatas@usc.edu \
    "job_scripts/sim_chr_rep_model.sh" \
        "${TREE_DIR}" \
        "${PICKLED_DEMO_META}" \
        "${VCF_DIR}" \
        "${PLINK_BED_DIR}" \
        "${POP_INFO_DIR}" \
        "${ANC_DIR}" \
        "${GLOBAL_ANC_DIR}" \
        "${SAMPLE_SIZE}" \
        "${NUM_REPS}" \
        "${MSPRIME_MODEL}" \
        "${GENETIC_MAP}" \
        "${GENERATION_TIME}" \
        "${MUTATION_RATE}" \
        "${T_AF_YEARS}" \
        "${T_OOA_YEARS}" \
        "${T_EU0_YEARS}" \
        "${T_EG_YEARS}" \
        "${R_EU0}" \
        "${R_EU}" \
        "${R_AF}" \
        "${N_A}" \
        "${N_AF1}" \
        "${N_B}" \
        "${N_EU0}" \
        "${M_AF_B}" \
        "${M_AF_EU}" \
        "${ADMIXTURE_TIME}" \
        "${ADMIX_GENERATION_COUNT}" \
        "${ADMIX_MIXING_GENERATION_COUNT}" \
        "${admix_ne_by_generation}" \
        "${admix_afr_props_by_generation}" \
        "${admix_eur_props_by_generation}" \
        "${admix_prioradmix_props_by_generation}" \
        "${ADMIX_MODERN_GROWTH_RATE}" \
        "${CENSUS_TIME_OFFSET}" \
        -- \
        "${CHROMS[@]}" \
        -- \
        "${POPS[@]}"
)
log_msg "submitted simulation/data-generation array; jid=${sim_jid}"

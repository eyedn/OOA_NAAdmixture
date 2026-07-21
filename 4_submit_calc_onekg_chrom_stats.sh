#!/usr/bin/env bash

###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           4_submit_calc_onekg_chrom_stats.sh
###############################################################################

# workflow: prepare and analyze 1000 Genomes chromosome-population inputs, then
# combine completed chromosome outputs after their dependent array jobs finish.

##### set up ##################################################################
set -euo pipefail

# determine repo location; all referenced scripts run from this checkout.
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# source shared empirical constants and timestamped logging.
source "${script_dir}/other_scripts/const.sh"
source "${script_dir}/other_scripts/log_msg.sh"


##### chromosome-population pipeline jobs #####################################
# create empirical output directories shared by preparation and statistics jobs.
mkdir -p "${ONEKG_OUT_VCF_DIR}" "${ONEKG_OUT_PLINK_BED_DIR}" \
    "${ONEKG_OUT_POP_INFO_DIR}" "${ONEKG_OUT_ADMIXTURE_DIR}" \
    "${ONEKG_OUT_KING_DIR}" "${ONEKG_OUT_STATS_DIR}"

# use chromosome-major, population-minor task IDs in both worker arrays.
array_size=$((${#CHROMS[@]} * ${#ONEKG_POPS[@]}))

# prepare one chromosome-population input subset per array task.
parse_jname="parse1kG"
mkdir -p "/home1/karatas/logs/${parse_jname}"
parse_jid=$(sbatch \
    --parsable \
    --chdir="${script_dir}" \
    --job-name="${parse_jname}" \
    --array="1-${array_size}%${MAX_JOBS}" \
    --cpus-per-task="${STATS_CPUS_PER_TASK}" \
    --mem="${STATS_MEM}" \
    --time=1-00:00:00 \
    --partition=qcb \
    --account=jazlynmo_738 \
    --nodes=1 \
    --output="/home1/karatas/logs/${parse_jname}/%A_%a.%x.out" \
    --error="/home1/karatas/logs/${parse_jname}/%A_%a.%x.err" \
    --mail-type=ALL \
    --mail-user=karatas@usc.edu \
    "job_scripts/parse_onekg_chr_pop.sh" \
        "${ONEKG_VCF_PREFIX}" \
        "${ONEKG_VCF_SUFFIX}" \
        "${ONEKG_UNRELS_FILE}" \
        "${ONEKG_PLINK_FAM_FILE}" \
        "${ONEKG_OUT_VCF_DIR}" \
        "${ONEKG_OUT_PLINK_BED_DIR}" \
        "${ONEKG_OUT_POP_INFO_DIR}" \
        -- \
        "${CHROMS[@]}" \
        -- \
        "${ONEKG_POPS[@]}"
)

# calculate statistics only after all preparation tasks complete successfully.
stats_jname="stats1kG"
mkdir -p "/home1/karatas/logs/${stats_jname}"
stats_jid=$(sbatch \
    --parsable \
    --dependency="afterok:${parse_jid}" \
    --chdir="${script_dir}" \
    --job-name="${stats_jname}" \
    --array="1-${array_size}%${MAX_JOBS}" \
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
    "job_scripts/calc_onekg_chr_pop_stats.sh" \
        "${ONEKG_VCF_PREFIX}" \
        "${ONEKG_VCF_SUFFIX}" \
        "${ONEKG_UNRELS_FILE}" \
        "${ONEKG_PLINK_FAM_FILE}" \
        "${ONEKG_CHR_LENS}" \
        "${ONEKG_OUT_VCF_DIR}" \
        "${ONEKG_OUT_PLINK_BED_DIR}" \
        "${ONEKG_OUT_POP_INFO_DIR}" \
        "${ONEKG_OUT_KING_DIR}" \
        "${ONEKG_OUT_STATS_DIR}" \
        "${MUTATION_RATE}" \
        -- \
        "${CHROMS[@]}" \
        -- \
        "${ONEKG_POPS[@]}"
)

##### chromosome statistics combination #######################################
# combine each chromosome only after its preparation and statistics pipeline.
comb_jname="comb1kG"
mkdir -p "/home1/karatas/logs/${comb_jname}"
comb_jid=$(sbatch \
    --parsable \
    --dependency="afterok:${stats_jid}" \
    --chdir="${script_dir}" \
    --job-name="${comb_jname}" \
    --array="1-${#CHROMS[@]}%${MAX_JOBS}" \
    --mem="${COMB_MEM}" \
    --time=1-00:00:00 \
    --partition=qcb \
    --account=jazlynmo_738 \
    --nodes=1 \
    --output="/home1/karatas/logs/${comb_jname}/%A_%a.%x.out" \
    --error="/home1/karatas/logs/${comb_jname}/%A_%a.%x.err" \
    --mail-type=ALL \
    --mail-user=karatas@usc.edu \
    "job_scripts/combine_onekg_chr_stats.sh" \
        "${ONEKG_UNRELS_FILE}" "${ONEKG_PLINK_FAM_FILE}" \
        "${ONEKG_OUT_VCF_DIR}" "${ONEKG_OUT_PLINK_BED_DIR}" \
        "${ONEKG_OUT_POP_INFO_DIR}" "${ONEKG_OUT_ADMIXTURE_DIR}" \
        "${ONEKG_OUT_STATS_DIR}" "${ADMIXTURE_LD_WINDOW}" \
        "${ADMIXTURE_LD_STEP}" "${ADMIXTURE_LD_R2}" \
        -- "${CHROMS[@]}" -- "${ONEKG_POPS[@]}"
)

log_msg "submitted 1000 Genomes chromosome pipeline; jid=${comb_jid}"

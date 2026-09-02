#!/usr/bin/env bash

###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           5_submit_calc_simOnekgDownsample_chrom_stats.sh
###############################################################################

# workflow: submit chromosome-statistics simulation jobs at 1kG SNP densities.


##### set up ##################################################################
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${script_dir}/other_scripts/const.sh"
source "${script_dir}/other_scripts/log_msg.sh"

if [[ ! "${OUT_DIR}" =~ _(small|large)$ ]]; then
    echo "ERROR: OUT_DIR must end in _small or _large: ${OUT_DIR}" >&2
    exit 1
fi
if [[ -z ${OOA_NAADMIXTURE_CONDA} \
    || -z ${FASTSTRUCTURE_CONDA_ENV} \
    || -z ${FASTSTRUCTURE_PRIOR} \
    || -z ${FASTSTRUCTURE_CV} \
    || ${#SIM_UNSUPERVISED_KS[@]} -eq 0 ]]; then
    echo "ERROR: inference environments and settings must be configured" >&2
    exit 1
fi
if [[ ! -x "${ADMIXTURE_EXEC}" ]]; then
    echo "ERROR: ADMIXTURE executable is not executable: ${ADMIXTURE_EXEC}" >&2
    exit 1
fi
for faststructure_script in \
    "${FASTSTRUCTURE_STRUCTURE_PY}" \
    "${FASTSTRUCTURE_CHOOSE_K_PY}"; do
    if [[ ! -f "${faststructure_script}" ]]; then
        echo "ERROR: missing fastStructure script ${faststructure_script}" >&2
        exit 1
    fi
done

mkdir -p \
    "${SIMONEKGDOWNSAMPLE_VCF_DIR}" \
    "${SIMONEKGDOWNSAMPLE_PLINK_BED_DIR}" \
    "${SIMONEKGDOWNSAMPLE_POP_INFO_DIR}" \
    "${SIMONEKGDOWNSAMPLE_ADMIXTURE_DIR}" \
    "${SIMONEKGDOWNSAMPLE_FASTSTRUCTURE_DIR}" \
    "${SIMONEKGDOWNSAMPLE_KING_DIR}" \
    "${SIMONEKGDOWNSAMPLE_STATS_DIR}"

array_size=$((NUM_REPS * ${#CHROMS[@]}))


##### density and input preparation ###########################################
prepare_jname="prepSim1kGDens"
mkdir -p "/home1/karatas/logs/${prepare_jname}"
prepare_jid=$(sbatch \
    --parsable \
    --chdir="${script_dir}" \
    --job-name="${prepare_jname}" \
    --array="1-${array_size}%${MAX_JOBS}" \
    --cpus-per-task="${SIMONEKGDOWNSAMPLE_CPUS_PER_TASK}" \
    --mem="${SIMONEKGDOWNSAMPLE_MEM}" \
    --time=1-00:00:00 \
    --partition="${PARTITION}" \
    --account="${ACCOUNT}" \
    --nodes=1 \
    --output="/home1/karatas/logs/${prepare_jname}/%A_%a.%x.out" \
    --error="/home1/karatas/logs/${prepare_jname}/%A_%a.%x.err" \
    --mail-type="${MAIL_TYPE}" \
    --mail-user="${MAIL_USER}" \
    "job_scripts/prepare_simOnekgDownsample_chr_rep.sh" \
        "${OOA_NAADMIXTURE_CONDA}" \
        "${VCF_DIR}" \
        "${POP_INFO_DIR}" \
        "${ONEKG_OUT_VCF_DIR}" \
        "${ONEKG_ADMIXED_POP_LABEL}" \
        "${SPAN_INCL_FILE}" \
        "${INTERGENTIC_FILE}" \
        "${SIMONEKGDOWNSAMPLE_VCF_DIR}" \
        "${SIMONEKGDOWNSAMPLE_PLINK_BED_DIR}" \
        "${SIMONEKGDOWNSAMPLE_POP_INFO_DIR}" \
        "${SIMONEKGDOWNSAMPLE_ADMIXTURE_DIR}" \
        "${SIMONEKGDOWNSAMPLE_FASTSTRUCTURE_DIR}" \
        "${SIMONEKGDOWNSAMPLE_KING_DIR}" \
        "${SIMONEKGDOWNSAMPLE_STATS_DIR}" \
        "${SNP_DENSITY_WINDOW_BP}" \
        "${NUM_REPS}" \
        "${GENETIC_MAP}" \
        -- \
        "${CHROMS[@]}" \
        -- \
        "${POPS[@]}")
log_msg "submitted simOnekgDownsample preparation; jid=${prepare_jid}"


##### chromosome statistics ##################################################
stats_jname="statsSim1kGDens"
mkdir -p "/home1/karatas/logs/${stats_jname}"
stats_jid=$(sbatch \
    --parsable \
    --dependency="afterok:${prepare_jid}" \
    --chdir="${script_dir}" \
    --job-name="${stats_jname}" \
    --array="1-${array_size}%${MAX_JOBS}" \
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
    "job_scripts/calc_simOnekgDownsample_chr_rep_stats.sh" \
        "${OOA_NAADMIXTURE_CONDA}" \
        "${ADMIXTURE_EXEC}" \
        "${FASTSTRUCTURE_CONDA_ENV}" \
        "${FASTSTRUCTURE_STRUCTURE_PY}" \
        "${FASTSTRUCTURE_CHOOSE_K_PY}" \
        "${FASTSTRUCTURE_PRIOR}" \
        "${FASTSTRUCTURE_CV}" \
        "${SIMONEKGDOWNSAMPLE_FASTSTRUCTURE_DIR}" \
        "${SIMONEKGDOWNSAMPLE_VCF_DIR}" \
        "${SIMONEKGDOWNSAMPLE_PLINK_BED_DIR}" \
        "${SIMONEKGDOWNSAMPLE_POP_INFO_DIR}" \
        "${SIMONEKGDOWNSAMPLE_ADMIXTURE_DIR}" \
        "${SIMONEKGDOWNSAMPLE_KING_DIR}" \
        "${SIMONEKGDOWNSAMPLE_STATS_DIR}" \
        "${ONEKG_CHR_LENS}" \
        "${INTERGENTIC_FILE}" \
        "${SPAN_INCL_FILE}" \
        "${NUM_REPS}" \
        "${GENETIC_MAP}" \
        "${MUTATION_RATE}" \
        "${KIN_CUTOFF}" \
        "${ADMIXTURE_LD_WINDOW}" \
        "${ADMIXTURE_LD_STEP}" \
        "${ADMIXTURE_LD_R2}" \
        "${LD_DECAY_WINDOW_SIZE_BP}" \
        "${LD_DECAY_DISTANCE_BIN_BP}" \
        "${LD_DECAY_MAF_THRESHOLD}" \
        "${SFS_SIZE}" \
        "${SFS_SIZE_POP_REF}" \
        -- \
        "${SIM_UNSUPERVISED_KS[@]}" \
        -- \
        "${CHROMS[@]}" \
        -- \
        "${POPS[@]}")
log_msg "submitted simOnekgDownsample statistics; jid=${stats_jid}"


##### chromosome statistics combination ######################################
combine_jname="combSim1kGDens"
mkdir -p "/home1/karatas/logs/${combine_jname}"
combine_jid=$(sbatch \
    --parsable \
    --dependency="afterok:${stats_jid}" \
    --chdir="${script_dir}" \
    --job-name="${combine_jname}" \
    --array="1-${#CHROMS[@]}%${MAX_JOBS}" \
    --cpus-per-task=1 \
    --mem="${ONEKG_COMB_MEM}" \
    --time=1-00:00:00 \
    --partition="${PARTITION}" \
    --account="${ACCOUNT}" \
    --nodes=1 \
    --output="/home1/karatas/logs/${combine_jname}/%A_%a.%x.out" \
    --error="/home1/karatas/logs/${combine_jname}/%A_%a.%x.err" \
    --mail-type="${MAIL_TYPE}" \
    --mail-user="${MAIL_USER}" \
    "job_scripts/combine_sim_chr_stats.sh" \
        "${OOA_NAADMIXTURE_CONDA}" \
        "${SIMONEKGDOWNSAMPLE_STATS_DIR}" \
        "${NUM_REPS}" \
        -- \
        "${CHROMS[@]}" \
        -- \
        ancestry_ADMIXTURE_super \
        ancestry_ADMIXTURE_multik \
        ancestry_fastStructure_multik \
        fastStructure_chooseK \
        kinship \
        kinship_unrelated \
        pi_theta_stats_intergenic \
        pi_theta_stats_full_callable_chrom \
        sfs \
        ld_decay \
        variant_qc \
        snp_density)
log_msg "submitted simOnekgDownsample combination; jid=${combine_jid}"

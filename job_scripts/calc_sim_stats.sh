#!/usr/bin/env bash

###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology 
#           Mooney Lab
#           ---
#           run_stats.sh
###############################################################################

#SBATCH --time=1-00:00:00
#SBATCH --partition=qcb
#SBATCH --account=jazlynmo_738
#SBATCH --nodes=1
#SBATCH --output=/home1/karatas/logs/calcOOANAA/calcOOANAA.%A_%a.%x.out
#SBATCH --error=/home1/karatas/logs/calcOOANAA/calcOOANAA.%A_%a.%x.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=karatas@usc.edu


set -euo pipefail

module purge
ml gcc/13.3.0 htslib/1.19.1 plink2/2.00a4.3 conda
source /apps/conda/miniforge3/25.3.0/etc/profile.d/conda.sh
conda activate OOA_NAAdmixture
export PATH="${HOME}/.conda/envs/OOA_NAAdmixture/bin:${PATH}"

project_dir="$(pwd)"
source "${project_dir}/other_scripts/log_msg.sh"

: "${SLURM_ARRAY_TASK_ID:?ERROR: run as a Slurm array}"


# input variables
tree_dir="$1"
plink_bed_dir="$2"
pop_info_dir="$3"
admixture_dir="$4"
global_anc_dir="$5"
king_dir="$6"
stats_dir="$7"
sample_size="$8"
num_reps="$9"
chr="${10}"
genetic_map="${11}"
mutation_rate="${12}"
kin_cutoff="${13}"
admixture_ld_window="${14}"
admixture_ld_step="${15}"
admixture_ld_r2="${16}"
shift 16
shift 1 # skip the "--" from input arguments
pops=( "$@" )

# derived variables
rep="${SLURM_ARRAY_TASK_ID}"
num_threads="${SLURM_CPUS_PER_TASK:-1}"
prefix="${genetic_map}_${rep}_all"
tree_tsz_path="${tree_dir}/${prefix}.ts.tsz"
plink_bed_prefix="${plink_bed_dir}/${prefix}"
sample_metadata_path="${pop_info_dir}/${genetic_map}_${rep}.sample_metadata.tsv"
admixture_prefix="${admixture_dir}/${prefix}"
unrelated_keep_path="${admixture_dir}/${prefix}.king_unrelated.keep"
ld_prune_prefix="${admixture_dir}/${prefix}.ld_prune"
if (( rep < 1 || rep > num_reps )); then
    echo "ERROR: invalid replicate ${rep}; expected 1..${num_reps}" >&2
    exit 1
fi

# calculated KING kinships coefficients
mkdir -p "${admixture_dir}" "${king_dir}" "${stats_dir}"
log_msg "running within-population KING for rep=${rep}"
: > "${unrelated_keep_path}"
for pop in "${pops[@]}"; do
    subset_path="${king_dir}/${genetic_map}_${rep}_${pop}.subset"
    out_prefix="${king_dir}/${genetic_map}_${rep}_${pop}"
    python "${project_dir}/job_scripts/write_pop_subset.py" \
        --subset-path "${subset_path}" \
        --pop "${pop}" \
        --sample-size "${sample_size}" \
        --pops "${pops[@]}"
    plink2 \
        --bfile "${plink_bed_prefix}" \
        --keep "${subset_path}" \
        --threads "${num_threads}" \
        --king-cutoff "${kin_cutoff}" \
        --make-king \
        --make-king-table \
        --out "${out_prefix}"
    retained_path="${out_prefix}.king.cutoff.in.id"
    if [[ ! -s "${retained_path}" ]]; then
        echo "ERROR: missing KING retained sample file ${retained_path}" >&2
        exit 1
    fi
    cat "${retained_path}" >> "${unrelated_keep_path}"
done

# create LD-pruned unrelated input for "admixture --supervised"
log_msg "LD-pruning unrelated ADMIXTURE samples for rep=${rep}"
plink2 \
    --bfile "${plink_bed_prefix}" \
    --keep "${unrelated_keep_path}" \
    --threads "${num_threads}" \
    --indep-pairwise "${admixture_ld_window}" \
        "${admixture_ld_step}" \
        "${admixture_ld_r2}" \
    --out "${ld_prune_prefix}"

if [[ ! -s "${ld_prune_prefix}.prune.in" ]]; then
    echo "ERROR: missing LD-pruned SNP list ${ld_prune_prefix}.prune.in" >&2
    exit 1
fi

log_msg "creating final ADMIXTURE BED for rep=${rep}"
plink2 \
    --bfile "${plink_bed_prefix}" \
    --keep "${unrelated_keep_path}" \
    --extract "${ld_prune_prefix}.prune.in" \
    --threads "${num_threads}" \
    --make-bed \
    --out "${admixture_prefix}"

if [[ ! -s "${admixture_prefix}.bed" || ! -s "${admixture_prefix}.bim" \
    || ! -s "${admixture_prefix}.fam" ]]; then
    echo "ERROR: failed to create final ADMIXTURE BED set" >&2
    exit 1
fi

log_msg "writing final supervised ADMIXTURE pop file for rep=${rep}"
python "${project_dir}/job_scripts/write_admixture_pop.py" \
    --sample-metadata-path "${sample_metadata_path}" \
    --fam-path "${admixture_prefix}.fam" \
    --pop-path "${admixture_prefix}.pop"

log_msg "running supervised ADMIXTURE for rep=${rep}"
(
    cd "${admixture_dir}"
    "${HOME}/software/ADMIXTURE/admixture_linux-1.4.0/admixture" \
        --supervised \
        -j"${num_threads}" \
        -s "${rep}" \
        "${prefix}.bed" \
        2
)

# calculate summaries on pi, theta, sfs, ld, and kinship
log_msg "parsing stats for rep=${rep}"
python "${project_dir}/job_scripts/calc_sim_stats.py" \
    --rep "${rep}" \
    --tree-tsz-path "${tree_tsz_path}" \
    --plink-bed-prefix "${plink_bed_prefix}" \
    --sample-metadata-path "${sample_metadata_path}" \
    --admixture-fam-path "${admixture_prefix}.fam" \
    --admixture-dir "${admixture_dir}" \
    --global-anc-dir "${global_anc_dir}" \
    --king-dir "${king_dir}" \
    --stats-dir "${stats_dir}" \
    --sample-size "${sample_size}" \
    --chr "${chr}" \
    --genetic-map "${genetic_map}" \
    --mutation-rate "${mutation_rate}" \
    --pops "${pops[@]}"

log_msg "done with statistics for rep=${rep}"

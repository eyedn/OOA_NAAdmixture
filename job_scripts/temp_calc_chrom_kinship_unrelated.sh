#!/usr/bin/env bash

###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           temp_calc_chrom_kinship_unrelated.sh
###############################################################################


set -euo pipefail

module purge
ml gcc/13.3.0 plink2/2.00a4.3 conda
source /apps/conda/miniforge3/25.3.0/etc/profile.d/conda.sh
conda activate OOA_NAAdmixture
export PATH="${HOME}/.conda/envs/OOA_NAAdmixture/bin:${PATH}"

project_dir="$(pwd)"
source "${project_dir}/other_scripts/log_msg.sh"
source "${project_dir}/other_scripts/map_slurm_task.sh"

: "${SLURM_ARRAY_TASK_ID:?ERROR: run as a Slurm array}"


# input variables
plink_bed_dir="$1"
king_dir="$2"
stats_dir="$3"
num_reps="$4"
genetic_map="$5"
shift 5
shift 1 # skip the "--" from input arguments
chroms=()
while [[ "$1" != "--" ]]; do
    chroms+=( "$1" )
    shift
done
shift 1 # skip the "--" from input arguments
pops=( "$@" )

# derived variables
read -r rep chr < <(
    map_slurm_task_to_rep_chr \
        "${SLURM_ARRAY_TASK_ID}" \
        "${num_reps}" \
        "${chroms[@]}"
)
num_threads="${SLURM_CPUS_PER_TASK:-1}"
plink_bed_prefix="${plink_bed_dir}/${genetic_map}_${rep}_chr${chr}_all"

# require existing PLINK inputs before running the focused backfill
for extension in bed bim fam; do
    bed_path="${plink_bed_prefix}.${extension}"
    if [[ ! -s "${bed_path}" ]]; then
        echo "ERROR: missing PLINK input ${bed_path}" >&2
        exit 1
    fi
done

mkdir -p "${king_dir}" "${stats_dir}"
log_msg "running unrelated-only chromosome KING rep=${rep} chr=${chr}"
for pop in "${pops[@]}"; do
    pop_prefix="${genetic_map}_${rep}_chr${chr}_${pop}"
    retained_path="${king_dir}/${pop_prefix}.king.cutoff.in.id"
    if [[ ! -s "${retained_path}" ]]; then
        echo "ERROR: missing KING retained sample file ${retained_path}" >&2
        exit 1
    fi
    unrelated_out_prefix="${king_dir}/${pop_prefix}_unrelated"
    plink2 \
        --bfile "${plink_bed_prefix}" \
        --keep "${retained_path}" \
        --threads "${num_threads}" \
        --make-king \
        --make-king-table \
        --out "${unrelated_out_prefix}"
done

python "${project_dir}/job_scripts/python_utils/write_unrelated_kinship.py" \
    --rep "${rep}" \
    --king-dir "${king_dir}" \
    --stats-dir "${stats_dir}" \
    --genetic-map "${genetic_map}" \
    --chr "${chr}" \
    --pops "${pops[@]}"

log_msg "done unrelated-only chromosome KING rep=${rep} chr=${chr}"

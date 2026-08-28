#!/usr/bin/env bash

###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           regenerate_onekg_chr_sfs.sh
###############################################################################

# workflow: rescan one empirical chromosome once and regenerate all SFS tables.


##### set up ##################################################################
set -euo pipefail
conda_env="$1"
vcf_prefix="$2"
vcf_suffix="$3"
unrels_path="$4"
fam_path="$5"
stats_dir="$6"
sfs_size="$7"
sfs_size_pop_ref="$8"
shift 8
shift
chroms=()
while [[ "$1" != "--" ]]; do
    chroms+=( "$1" )
    shift
done
shift
pops=( "$@" )

module purge
ml gcc/13.3.0 htslib/1.19.1 conda
source /apps/conda/miniforge3/25.3.0/etc/profile.d/conda.sh
conda activate "${conda_env}"
project_dir="$(pwd)"
source "${project_dir}/other_scripts/log_msg.sh"
: "${SLURM_ARRAY_TASK_ID:?ERROR: run as a Slurm array}"
chrom_index=$((SLURM_ARRAY_TASK_ID - 1))
if (( chrom_index < 0 || chrom_index >= ${#chroms[@]} )); then
    echo "ERROR: invalid chromosome task ID ${SLURM_ARRAY_TASK_ID}" >&2
    exit 1
fi
chrom="${chroms[chrom_index]}"


##### regenerate SFS ##########################################################
log_msg "regenerating 1000 Genomes SFS chr=${chrom}"
python "${project_dir}/job_scripts/python_utils/regenerate_onekg_chr_sfs.py" \
    --vcf-path "${vcf_prefix}${chrom}${vcf_suffix}" \
    --unrels-path "${unrels_path}" \
    --fam-path "${fam_path}" \
    --stats-dir "${stats_dir}" \
    --chrom "${chrom}" \
    --sfs-size "${sfs_size}" \
    --sfs-size-pop-ref "${sfs_size_pop_ref}" \
    --pops "${pops[@]}"
log_msg "completed 1000 Genomes SFS chr=${chrom}"

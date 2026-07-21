#!/usr/bin/env bash

###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           calc_onekg_genome_pop_stats.sh
###############################################################################

# workflow: merge 1000 Genomes chromosomes and calculate genome statistics.


##### set up ##################################################################
set -euo pipefail

# load the software environment used for PLINK and Python statistics.
module purge
ml gcc/13.3.0 plink2/2.00a4.3 conda
source /apps/conda/miniforge3/25.3.0/etc/profile.d/conda.sh
conda activate OOA_NAAdmixture

# load shared logging and require a population-indexed Slurm task.
project_dir="$(pwd)"
source "${project_dir}/other_scripts/log_msg.sh"
: "${SLURM_ARRAY_TASK_ID:?ERROR: run as a Slurm array}"


##### variables ###############################################################
# read fixed inputs followed by chromosome and population arrays.
bed_dir="$1"
pop_info_dir="$2"
king_dir="$3"
stats_dir="$4"
shift 4
shift
chroms=()
while [[ "$1" != "--" ]]; do
    chroms+=( "$1" )
    shift
done
shift
pops=( "$@" )

# map the one-dimensional task array directly to the population list.
pop_index=$((SLURM_ARRAY_TASK_ID - 1))
if (( pop_index < 0 || pop_index >= ${#pops[@]} )); then
    echo "ERROR: invalid population task ID ${SLURM_ARRAY_TASK_ID}" >&2
    exit 1
fi
pop="${pops[pop_index]}"
num_threads="${SLURM_CPUS_PER_TASK:-1}"
merge_list="${bed_dir}/onekg.rep_0.genome.${pop}.merge_list"
genome_prefix="${bed_dir}/onekg.rep_0.genome.${pop}"
king_prefix="${king_dir}/onekg.rep_0.genome.${pop}"
pop_keep="${pop_info_dir}/onekg.rep_0.chr${chroms[0]}.${pop}.keep"


##### genome PLINK and KING outputs ###########################################
# merge chromosome PLINK sets, then calculate KING for the target population.
: > "${merge_list}"
for chr in "${chroms[@]}"; do
    printf '%s\n' "${bed_dir}/onekg.rep_0.chr${chr}.all" \
        >> "${merge_list}"
done
plink2 \
    --pmerge-list "${merge_list}" bfile \
    --threads "${num_threads}" \
    --make-bed \
    --out "${genome_prefix}"
plink2 \
    --bfile "${genome_prefix}" \
    --keep "${pop_keep}" \
    --threads "${num_threads}" \
    --make-king-table \
    --out "${king_prefix}"


##### statistics ##############################################################
# aggregate chromosome statistics using the genome-level KING output.
python "${project_dir}/job_scripts/python_utils/calc_onekg_stats.py" \
    --analysis-level "genome" \
    --king-path "${king_prefix}.kin0" \
    --stats-dir "${stats_dir}" \
    --pop "${pop}" \
    --chroms "${chroms[@]}"
log_msg "done with 1000 Genomes genome statistics pop=${pop}"

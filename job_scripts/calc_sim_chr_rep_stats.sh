#!/usr/bin/env bash

###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           calc_sim_chr_rep_stats.sh
###############################################################################

# workflow: calculate simulation statistics for one chromosome and replicate.


##### set up ##################################################################
set -euo pipefail

# load required HPC modules and conda env
module purge
ml gcc/13.3.0 htslib/1.19.1 plink2/2.00a4.3 conda
source /apps/conda/miniforge3/25.3.0/etc/profile.d/conda.sh
conda activate OOA_NAAdmixture
export PATH="${HOME}/.conda/envs/OOA_NAAdmixture/bin:${PATH}"

# load shared functions; note, all scripts should exist in the execution repo
project_dir="$(pwd)"
source "${project_dir}/other_scripts/log_msg.sh"
source "${project_dir}/other_scripts/map_slurm_task.sh"

# require execution as a Slurm array task
: "${SLURM_ARRAY_TASK_ID:?ERROR: run as a Slurm array}"


##### variables ###############################################################
# read fixed statistics inputs followed by chromosome and population arrays.
tree_dir="$1"
plink_bed_dir="$2"
pop_info_dir="$3"
admixture_dir="$4"
global_anc_dir="$5"
king_dir="$6"
stats_dir="$7"
sample_size="$8"
num_reps="$9"
genetic_map="${10}"
mutation_rate="${11}"
kin_cutoff="${12}"
admixture_ld_window="${13}"
admixture_ld_step="${14}"
admixture_ld_r2="${15}"
shift 15
shift 1 # skip the "--" from input arguments
chroms=()
while [[ "$1" != "--" ]]; do
    chroms+=( "$1" )
    shift
done
shift 1 # skip the "--" from input arguments
pops=( "$@" )

# map the array-task ID to a chrom and rep; note, tasks are ordered
# chromosome-major, with replicate-minor indexing
read -r chr rep < <(
    map_slurm_task_to_chr_rep \
        "${SLURM_ARRAY_TASK_ID}" \
        "${num_reps}" \
        -- \
        "${chroms[@]}"
)

# derive task-specific input, output, and intermediate paths.
num_threads="${SLURM_CPUS_PER_TASK:-1}"
prefix="${genetic_map}_${rep}_chr${chr}_all"
tree_tsz_path="${tree_dir}/${prefix}.ts.tsz"
plink_bed_prefix="${plink_bed_dir}/${prefix}"
sample_metadata_path="${pop_info_dir}/${genetic_map}_${rep}_chr${chr}.sample_metadata.tsv"
admixture_prefix="${admixture_dir}/${prefix}"
unrelated_keep_path="${admixture_dir}/${prefix}.king_unrelated.keep"
ld_prune_prefix="${admixture_dir}/${prefix}.ld_prune"


##### KING coefficients #######################################################
# create all required output directories
mkdir -p "${admixture_dir}" "${king_dir}" "${stats_dir}"

# calculate relatedness and unrelated-sample handoffs for each population.
log_msg "running within-population KING for rep=${rep} chr=${chr}"
: > "${unrelated_keep_path}"
for pop in "${pops[@]}"; do
    subset_path="${king_dir}/${genetic_map}_${rep}_chr${chr}_${pop}.subset"
    out_prefix="${king_dir}/${genetic_map}_${rep}_chr${chr}_${pop}"

    # create the population-specific sample subset for PLINK.
    python "${project_dir}/job_scripts/python_utils/write_pop_subset.py" \
        --subset-path "${subset_path}" \
        --pop "${pop}" \
        --sample-size "${sample_size}" \
        --pops "${pops[@]}"

    # use population-specific sample IDs to calculate KING coefficients.
    plink2 \
        --bfile "${plink_bed_prefix}" \
        --keep "${subset_path}" \
        --threads "${num_threads}" \
        --king-cutoff "${kin_cutoff}" \
        --make-king \
        --make-king-table \
        --out "${out_prefix}"
    retained_path="${out_prefix}.king.cutoff.in.id"

    # verify KING output containing IDs retained by the kinship cutoff.
    if [[ ! -s "${retained_path}" ]]; then
        echo "ERROR: missing KING retained sample file ${retained_path}" >&2
        exit 1
    fi

    # rerun KING on unrelated subsets and combine retained IDs into one
    # single unrelated-sample keep file
    unrelated_out_prefix="${out_prefix}_unrelated"
    plink2 \
        --bfile "${plink_bed_prefix}" \
        --keep "${retained_path}" \
        --threads "${num_threads}" \
        --make-king \
        --make-king-table \
        --out "${unrelated_out_prefix}"
    cat "${retained_path}" >> "${unrelated_keep_path}"
done

# write unrelated-sample kinship summaries for the chromosome and replicate
# with all populations included
python "${project_dir}/job_scripts/python_utils/write_unrelated_kinship.py" \
    --rep "${rep}" \
    --king-dir "${king_dir}" \
    --stats-dir "${stats_dir}" \
    --genetic-map "${genetic_map}" \
    --chr "${chr}" \
    --pops "${pops[@]}"


##### ADMIXTURE ###############################################################
# LD-prune unrelated samples using the ADMIXTURE-recommended PLINK window,
# step size, and correlation threshold inputs.
log_msg "LD-pruning unrelated ADMIXTURE samples for rep=${rep} chr=${chr}"
plink2 \
    --bfile "${plink_bed_prefix}" \
    --keep "${unrelated_keep_path}" \
    --threads "${num_threads}" \
    --indep-pairwise "${admixture_ld_window}" \
        "${admixture_ld_step}" \
        "${admixture_ld_r2}" \
    --out "${ld_prune_prefix}"

# verify that the LD-pruned SNP list was created.
if [[ ! -s "${ld_prune_prefix}.prune.in" ]]; then
    echo "ERROR: missing LD-pruned SNP list ${ld_prune_prefix}.prune.in" >&2
    exit 1
fi

# create ADMIXTURE PLINK binaries from unrelated samples and LD-pruned SNPs.
log_msg "creating final ADMIXTURE BED for rep=${rep} chr=${chr}"
plink2 \
    --bfile "${plink_bed_prefix}" \
    --keep "${unrelated_keep_path}" \
    --extract "${ld_prune_prefix}.prune.in" \
    --threads "${num_threads}" \
    --make-bed \
    --out "${admixture_prefix}"

# verify that all ADMIXTURE PLINK binaries were created.
if [[ ! -s "${admixture_prefix}.bed" || ! -s "${admixture_prefix}.bim" \
    || ! -s "${admixture_prefix}.fam" ]]; then
    echo "ERROR: failed to create final ADMIXTURE BED set" >&2
    exit 1
fi

# generate the supervised ADMIXTURE .pop file from simulation metadata and
# final .fam order.
log_msg "writing final supervised ADMIXTURE pop file for rep=${rep} chr=${chr}"
python "${project_dir}/job_scripts/python_utils/write_admixture_pop.py" \
    --sample-metadata-path "${sample_metadata_path}" \
    --fam-path "${admixture_prefix}.fam" \
    --pop-path "${admixture_prefix}.pop"

# run supervised ADMIXTURE.
log_msg "running supervised ADMIXTURE for rep=${rep} chr=${chr}"
(
    cd "${admixture_dir}"
    "${HOME}/software/ADMIXTURE/admixture_linux-1.4.0/admixture" \
        --supervised \
        -j"${num_threads}" \
        -s "${rep}" \
        "${prefix}.bed" \
        2
)


##### statistics ##############################################################
# calculate chromosome ancestry, kinship, pi, theta, SFS, and LD summaries.
log_msg "parsing stats for rep=${rep} chr=${chr}"
python "${project_dir}/job_scripts/python_utils/calc_sim_stats.py" \
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

log_msg "done with statistics for rep=${rep} chr=${chr}"

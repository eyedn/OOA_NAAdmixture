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

conda_env="$1"
admixture_exec="$2"
faststructure_conda_env="$3"
faststructure_structure_py="$4"
faststructure_choose_k_py="$5"
faststructure_prior="$6"
faststructure_cv="$7"
faststructure_dir="$8"
shift 8

# load required HPC modules and conda env
module purge
ml gcc/13.3.0 htslib/1.19.1 plink2/2.00a4.3 conda
source /apps/conda/miniforge3/25.3.0/etc/profile.d/conda.sh
conda activate "${conda_env}"

# load shared functions; note, all scripts should exist in the execution repo
project_dir="$(pwd)"
source "${project_dir}/other_scripts/log_msg.sh"
source "${project_dir}/other_scripts/map_slurm_task.sh"

# require execution as a Slurm array task
: "${SLURM_ARRAY_TASK_ID:?ERROR: run as a Slurm array}"
if [[ ! -x "${admixture_exec}" ]]; then
    echo "ERROR: ADMIXTURE executable is not executable: ${admixture_exec}" >&2
    exit 1
fi
for faststructure_script in \
    "${faststructure_structure_py}" \
    "${faststructure_choose_k_py}"; do
    if [[ ! -f "${faststructure_script}" ]]; then
        echo "ERROR: missing fastStructure script ${faststructure_script}" >&2
        exit 1
    fi
done


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
ld_decay_window_size_bp="${16}"
ld_decay_distance_bin_bp="${17}"
ld_decay_maf_threshold="${18}"
sfs_size="${19}"
shift 19
shift 1 # skip the "--" from input arguments
unsupervised_ks=()
while [[ "$1" != "--" ]]; do
    unsupervised_ks+=( "$1" )
    shift
done
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
sample_metadata_path="${pop_info_dir}/${genetic_map}_${rep}"
sample_metadata_path+="_chr${chr}.sample_metadata.tsv"
admixture_prefix="${admixture_dir}/${prefix}"
faststructure_prefix="${faststructure_dir}/${prefix}"
faststructure_choose_k_path="${faststructure_prefix}.chooseK.txt"
unrelated_keep_path="${admixture_dir}/${prefix}.king_unrelated.keep"
ld_prune_prefix="${admixture_dir}/${prefix}.ld_prune"
inference_seed="$((1000 * rep + chr))"


##### KING coefficients #######################################################
# create all required output directories
mkdir -p "${admixture_dir}" "${faststructure_dir}" \
    "${king_dir}" "${stats_dir}"

# calculate relatedness and unrelated-sample handoffs for each population.
log_msg "running within-population KING for rep=${rep} chr=${chr}"
: > "${unrelated_keep_path}"
for pop in "${pops[@]}"; do
    subset_path="${king_dir}/${genetic_map}_${rep}_chr${chr}_${pop}.subset"
    out_prefix="${king_dir}/${genetic_map}_${rep}_chr${chr}_${pop}"

    # create the population-specific sample subset for PLINK.
    python "${project_dir}/job_scripts/python_utils/write_sim_pop_subset.py" \
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
python \
    "${project_dir}/job_scripts/python_utils/write_sim_unrelated_kinship.py" \
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
python "${project_dir}/job_scripts/python_utils/write_sim_admixture_pop.py" \
    --sample-metadata-path "${sample_metadata_path}" \
    --fam-path "${admixture_prefix}.fam" \
    --pop-path "${admixture_prefix}.pop"

# run supervised ADMIXTURE, preserve its K=2 Q file, then run each
# unsupervised K without allowing the K=2 artifact names to collide.
log_msg "running ADMIXTURE for rep=${rep} chr=${chr}"
(
    cd "${admixture_dir}"
    "${admixture_exec}" \
        --supervised \
        -j"${num_threads}" \
        -s "${rep}" \
        "${prefix}.bed" \
        2
    cp "${prefix}.2.Q" "${prefix}.supervised.2.Q"
    for k in "${unsupervised_ks[@]}"; do
        "${admixture_exec}" \
            -j"${num_threads}" \
            -s "${inference_seed}" \
            "${prefix}.bed" \
            "${k}"
        cp "${prefix}.${k}.Q" \
            "${prefix}.unsupervised.${k}.Q"
    done
)
if [[ ! -s "${admixture_prefix}.supervised.2.Q" ]]; then
    echo "ERROR: missing supervised ADMIXTURE Q file" >&2
    exit 1
fi
admixture_q_args=()
for k in "${unsupervised_ks[@]}"; do
    q_path="${admixture_prefix}.unsupervised.${k}.Q"
    if [[ ! -s "${q_path}" ]]; then
        echo "ERROR: missing unsupervised ADMIXTURE Q file ${q_path}" >&2
        exit 1
    fi
    admixture_q_args+=("--admixture-q-path" "${k}=${q_path}")
done

# run only fastStructure's Python 2 programs in its dedicated environment.
log_msg "deactivate ${conda_env}; activate ${faststructure_conda_env}"
conda deactivate
conda activate "${faststructure_conda_env}"
(
    for k in "${unsupervised_ks[@]}"; do
        python "${faststructure_structure_py}" \
            -K "${k}" \
            --input="${admixture_prefix}" \
            --output="${faststructure_prefix}" \
            --prior="${faststructure_prior}" \
            --cv="${faststructure_cv}" \
            --seed="${inference_seed}"
    done
    python "${faststructure_choose_k_py}" \
        --input="${faststructure_prefix}" \
        | tee "${faststructure_choose_k_path}"
)
log_msg "deactivate ${faststructure_conda_env}; activate ${conda_env}"
conda deactivate
conda activate "${conda_env}"

# validate fastStructure artifacts after restoring the Python 3 environment.
faststructure_q_args=()
for k in "${unsupervised_ks[@]}"; do
    q_path="${faststructure_prefix}.${k}.meanQ"
    for suffix in meanQ meanP log; do
        output_path="${faststructure_prefix}.${k}.${suffix}"
        if [[ ! -s "${output_path}" ]]; then
            echo "ERROR: missing fastStructure output ${output_path}" >&2
            exit 1
        fi
    done
    faststructure_q_args+=("--faststructure-q-path" "${k}=${q_path}")
done
if [[ ! -s "${faststructure_choose_k_path}" ]]; then
    echo "ERROR: missing fastStructure chooseK report" >&2
    exit 1
fi

# write neutral multi-K tables in final FAM order under Python 3.
python \
    "${project_dir}/job_scripts/python_utils/build_sim_inference_ancestry.py" \
    --rep "${rep}" \
    --sample-metadata-path "${sample_metadata_path}" \
    --admixture-fam-path "${admixture_prefix}.fam" \
    --supervised-q-path "${admixture_prefix}.supervised.2.Q" \
    "${admixture_q_args[@]}" \
    "${faststructure_q_args[@]}" \
    --faststructure-choose-k-path "${faststructure_choose_k_path}" \
    --faststructure-prior "${faststructure_prior}" \
    --faststructure-seed "${inference_seed}" \
    --stats-dir "${stats_dir}" \
    --chrom "${chr}"


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
    --sfs-size "${sfs_size}" \
    --chr "${chr}" \
    --genetic-map "${genetic_map}" \
    --mutation-rate "${mutation_rate}" \
    --ld-decay-window-size-bp "${ld_decay_window_size_bp}" \
    --ld-decay-distance-bin-bp "${ld_decay_distance_bin_bp}" \
    --ld-decay-maf-threshold "${ld_decay_maf_threshold}" \
    --pops "${pops[@]}"

log_msg "done with statistics for rep=${rep} chr=${chr}"

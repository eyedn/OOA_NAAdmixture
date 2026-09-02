#!/usr/bin/env bash

###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           calc_simOnekgDownsample_chr_rep_stats.sh
###############################################################################

# workflow: calculate empirical-style statistics for one density-matched task.


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

module purge
ml gcc/13.3.0 htslib/1.19.1 plink2/2.00a4.3 conda
source /apps/conda/miniforge3/25.3.0/etc/profile.d/conda.sh
conda activate "${conda_env}"

project_dir="$(pwd)"
source "${project_dir}/other_scripts/log_msg.sh"
source "${project_dir}/other_scripts/map_slurm_task.sh"
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
# read fixed statistics inputs followed by K, chromosome, and population arrays.
simOnekgDownsample_vcf_dir="$1"
simOnekgDownsample_plink_bed_dir="$2"
simOnekgDownsample_pop_info_dir="$3"
simOnekgDownsample_admixture_dir="$4"
simOnekgDownsample_king_dir="$5"
simOnekgDownsample_stats_dir="$6"
chr_lens_path="$7"
intergenic_file="$8"
span_incl_file="$9"
num_reps="${10}"
genetic_map="${11}"
mutation_rate="${12}"
kin_cutoff="${13}"
admixture_ld_window="${14}"
admixture_ld_step="${15}"
admixture_ld_r2="${16}"
ld_decay_window_size_bp="${17}"
ld_decay_distance_bin_bp="${18}"
ld_decay_maf_threshold="${19}"
sfs_size="${20}"
sfs_size_pop_ref="${21}"
shift 21
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
metadata_prefix="${genetic_map}_${rep}_chr${chr}"
simOnekgDownsample_vcf="${simOnekgDownsample_vcf_dir}/${prefix}.biallelic_snps.vcf.gz"
sim_metadata="${simOnekgDownsample_pop_info_dir}/${metadata_prefix}.sample_metadata.tsv"
all_bed_prefix="${simOnekgDownsample_plink_bed_dir}/${prefix}"
admixture_prefix="${simOnekgDownsample_admixture_dir}/${prefix}"
prune_prefix="${simOnekgDownsample_admixture_dir}/${prefix}.ld_prune"
faststructure_prefix="${faststructure_dir}/${prefix}"
faststructure_choose_k_path="${faststructure_prefix}.chooseK.txt"
unrelated_keep_path="${simOnekgDownsample_admixture_dir}/${prefix}.king_unrelated.keep"
inference_seed=$((1000 * rep + chr))

mkdir -p \
    "${simOnekgDownsample_admixture_dir}" \
    "${faststructure_dir}" \
    "${simOnekgDownsample_king_dir}" \
    "${simOnekgDownsample_stats_dir}"

for required_path in \
    "${simOnekgDownsample_vcf}" \
    "${sim_metadata}" \
    "${all_bed_prefix}.bed" \
    "${all_bed_prefix}.bim" \
    "${all_bed_prefix}.fam" \
    "${chr_lens_path}" \
    "${intergenic_file}" \
    "${span_incl_file}"; do
    if [[ ! -s "${required_path}" ]]; then
        echo "ERROR: missing required input ${required_path}" >&2
        exit 1
    fi
done


##### population statistics ##################################################
: > "${unrelated_keep_path}"
for pop in "${pops[@]}"; do
    pop_prefix="${metadata_prefix}_${pop}"
    pop_bed_prefix="${simOnekgDownsample_plink_bed_dir}/${pop_prefix}"
    intergenic_vcf="${simOnekgDownsample_vcf_dir}/${pop_prefix}.intergenic.vcf.gz"
    king_prefix="${simOnekgDownsample_king_dir}/${pop_prefix}"

    log_msg "running KING simOnekgDownsample rep=${rep} chr=${chr} pop=${pop}"
    plink2 \
        --bfile "${pop_bed_prefix}" \
        --threads "${num_threads}" \
        --king-cutoff "${kin_cutoff}" \
        --make-king \
        --make-king-table \
        --out "${king_prefix}"
    if [[ ! -s "${king_prefix}.kin0" ]]; then
        echo "ERROR: missing KING table ${king_prefix}.kin0" >&2
        exit 1
    fi
    retained_path="${king_prefix}.king.cutoff.in.id"
    if [[ ! -s "${retained_path}" ]]; then
        echo "ERROR: missing KING retained IDs ${retained_path}" >&2
        exit 1
    fi
    unrelated_king_prefix="${king_prefix}_unrelated"
    plink2 \
        --bfile "${pop_bed_prefix}" \
        --keep "${retained_path}" \
        --threads "${num_threads}" \
        --make-king \
        --make-king-table \
        --out "${unrelated_king_prefix}"
    if [[ ! -s "${unrelated_king_prefix}.kin0" ]]; then
        echo "ERROR: missing unrelated KING table" >&2
        exit 1
    fi
    cat "${retained_path}" >> "${unrelated_keep_path}"

    python "${project_dir}/job_scripts/python_utils/calc_onekg_stats.py" \
        --analysis-level chromosome \
        --rep "${rep}" \
        --vcf-path "${simOnekgDownsample_vcf}" \
        --ld-vcf-path "${simOnekgDownsample_vcf}" \
        --intergenic-vcf-path "${intergenic_vcf}" \
        --intergenic-bed-path "${intergenic_file}" \
        --span-incl-bed-path "${span_incl_file}" \
        --unrels-path "${sim_metadata}" \
        --fam-path "${all_bed_prefix}.fam" \
        --chr-lens-path "${chr_lens_path}" \
        --king-path "${unrelated_king_prefix}.kin0" \
        --stats-dir "${simOnekgDownsample_stats_dir}" \
        --chrom "${chr}" \
        --pop "${pop}" \
        --mutation-rate "${mutation_rate}" \
        --ld-decay-window-size-bp "${ld_decay_window_size_bp}" \
        --ld-decay-distance-bin-bp "${ld_decay_distance_bin_bp}" \
        --ld-decay-maf-threshold "${ld_decay_maf_threshold}" \
        --sfs-size "${sfs_size}" \
        --sfs-size-pop-ref "${sfs_size_pop_ref}" \
        --pops "${pops[@]}"
done

python \
"${project_dir}"/job_scripts/python_utils/write_simOnekgDownsample_kinship.py \
    --rep "${rep}" \
    --king-dir "${simOnekgDownsample_king_dir}" \
    --stats-dir "${simOnekgDownsample_stats_dir}" \
    --genetic-map "${genetic_map}" \
    --chrom "${chr}" \
    --pops "${pops[@]}"


##### ancestry inference #####################################################
log_msg "LD-pruning simOnekgDownsample rep=${rep} chr=${chr}"
plink2 \
    --bfile "${all_bed_prefix}" \
    --keep "${unrelated_keep_path}" \
    --indep-pairwise \
        "${admixture_ld_window}" \
        "${admixture_ld_step}" \
        "${admixture_ld_r2}" \
    --threads "${num_threads}" \
    --out "${prune_prefix}"
if [[ ! -s "${prune_prefix}.prune.in" ]]; then
    echo "ERROR: missing LD-pruned SNP list ${prune_prefix}.prune.in" >&2
    exit 1
fi
plink2 \
    --bfile "${all_bed_prefix}" \
    --keep "${unrelated_keep_path}" \
    --extract "${prune_prefix}.prune.in" \
    --threads "${num_threads}" \
    --make-bed \
    --out "${admixture_prefix}"
for extension in bed bim fam; do
    if [[ ! -s "${admixture_prefix}.${extension}" ]]; then
        echo "ERROR: missing ADMIXTURE PLINK ${extension} output" >&2
        exit 1
    fi
done

python "${project_dir}/job_scripts/python_utils/write_sim_admixture_pop.py" \
    --sample-metadata-path "${sim_metadata}" \
    --fam-path "${admixture_prefix}.fam" \
    --pop-path "${admixture_prefix}.pop"

admixture_q_args=()
(
    cd "${simOnekgDownsample_admixture_dir}"
    "${admixture_exec}" \
        --supervised \
        -j"${num_threads}" \
        -s "${inference_seed}" \
        "${prefix}.bed" \
        2
    cp "${prefix}.2.Q" "${prefix}.supervised.2.Q"
    for k in "${unsupervised_ks[@]}"; do
        "${admixture_exec}" \
            -j"${num_threads}" \
            -s "${inference_seed}" \
            "${prefix}.bed" \
            "${k}"
        cp "${prefix}.${k}.Q" "${prefix}.unsupervised.${k}.Q"
    done
)
if [[ ! -s "${admixture_prefix}.supervised.2.Q" ]]; then
    echo "ERROR: missing supervised ADMIXTURE Q output" >&2
    exit 1
fi
for k in "${unsupervised_ks[@]}"; do
    q_path="${admixture_prefix}.unsupervised.${k}.Q"
    if [[ ! -s "${q_path}" ]]; then
        echo "ERROR: missing ADMIXTURE output ${q_path}" >&2
        exit 1
    fi
    admixture_q_args+=( --admixture-q-path "${k}=${q_path}" )
done

log_msg "running fastStructure simOnekgDownsample rep=${rep} chr=${chr}"
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
conda deactivate
conda activate "${conda_env}"

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
    faststructure_q_args+=( --faststructure-q-path "${k}=${q_path}" )
done
if [[ ! -s "${faststructure_choose_k_path}" ]]; then
    echo "ERROR: missing fastStructure choose-K output" >&2
    exit 1
fi

python \
    "${project_dir}/job_scripts/python_utils/build_sim_inference_ancestry.py" \
    --rep "${rep}" \
    --sample-metadata-path "${sim_metadata}" \
    --admixture-fam-path "${admixture_prefix}.fam" \
    --supervised-q-path "${admixture_prefix}.supervised.2.Q" \
    "${admixture_q_args[@]}" \
    "${faststructure_q_args[@]}" \
    --faststructure-choose-k-path "${faststructure_choose_k_path}" \
    --faststructure-prior "${faststructure_prior}" \
    --faststructure-seed "${inference_seed}" \
    --stats-dir "${simOnekgDownsample_stats_dir}" \
    --chrom "${chr}"


##### replicate combination #################################################
python \
    "${project_dir}/job_scripts/python_utils/combine_onekg_chr_stats.py" \
    --stats-dir "${simOnekgDownsample_stats_dir}" \
    --chrom "${chr}" \
    --rep "${rep}" \
    --replicate-only \
    --preserve-kinship \
    --pops "${pops[@]}"

log_msg "done with simOnekgDownsample stats rep=${rep} chr=${chr}"

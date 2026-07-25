#!/usr/bin/env bash

###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           combine_onekg_genome_stats.sh
###############################################################################

# workflow: combine 1000 Genomes genome statistics into final table outputs.


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

# load the software environment used for PLINK, ADMIXTURE, and Python.
module purge
ml gcc/13.3.0 plink2/2.00a4.3 conda
source /apps/conda/miniforge3/25.3.0/etc/profile.d/conda.sh
conda activate "${conda_env}"

# load shared timestamped logging.
project_dir="$(pwd)"
source "${project_dir}/other_scripts/log_msg.sh"


##### variables ###############################################################
# read fixed inputs followed by K, chromosome, and population arrays.
unrels_path="$1"
source_fam_path="$2"
bed_dir="$3"
admixture_dir="$4"
stats_dir="$5"
ld_window="$6"
ld_step="$7"
ld_r2="$8"
shift 8
shift
unsupervised_ks=()
while [[ "$1" != "--" ]]; do
    unsupervised_ks+=( "$1" )
    shift
done
if (( ${#unsupervised_ks[@]} == 0 )); then
    echo "ERROR: at least one unsupervised K value is required" >&2
    exit 1
fi
shift
chroms=()
while [[ "$1" != "--" ]]; do
    chroms+=( "$1" )
    shift
done
shift
pops=( "$@" )

# derive genome-wide merge, pruning, and ADMIXTURE prefixes.
num_threads="${SLURM_CPUS_PER_TASK:-1}"
merge_list="${bed_dir}/onekg.rep_0.genome.merge_list"
genome_prefix="${bed_dir}/onekg.rep_0.genome.all"
prune_prefix="${admixture_dir}/onekg.rep_0.genome.ld_prune"
admixture_prefix="${admixture_dir}/onekg.rep_0.genome.all"
faststructure_prefix="${faststructure_dir}/onekg.rep_0.genome.all"
faststructure_choose_k_path="${faststructure_prefix}.chooseK.txt"
faststructure_seed="0"


##### global ancestry inference ###############################################
# merge and LD-prune the genome-wide input reused by both inference tools.
mkdir -p "${admixture_dir}" "${faststructure_dir}" "${stats_dir}"
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
    --indep-pairwise "${ld_window}" "${ld_step}" "${ld_r2}" \
    --threads "${num_threads}" \
    --out "${prune_prefix}"
if [[ ! -s "${prune_prefix}.prune.in" ]]; then
    echo "ERROR: missing genome LD-pruned SNP list" >&2
    exit 1
fi
plink2 \
    --bfile "${genome_prefix}" \
    --extract "${prune_prefix}.prune.in" \
    --threads "${num_threads}" \
    --make-bed \
    --out "${admixture_prefix}"
if [[ ! -s "${admixture_prefix}.bed" \
    || ! -s "${admixture_prefix}.bim" \
    || ! -s "${admixture_prefix}.fam" ]]; then
    echo "ERROR: failed to create final genome ADMIXTURE BED set" >&2
    exit 1
fi
# label reference samples and leave the admixed population unsupervised.
python "${project_dir}/job_scripts/python_utils/write_onekg_admixture_pop.py" \
    --unrels-path "${unrels_path}" \
    --source-fam-path "${source_fam_path}" \
    --admixture-fam-path "${admixture_prefix}.fam" \
    --pop-path "${admixture_prefix}.pop" \
    --afr-pop "${pops[0]}" \
    --eur-pop "${pops[1]}" \
    --admixed-pop "${pops[2]}"
(
    cd "${admixture_dir}"
    "${admixture_exec}" \
        --supervised \
        -j"${num_threads}" \
        "onekg.rep_0.genome.all.bed" 2
    cp "onekg.rep_0.genome.all.2.Q" \
        "onekg.rep_0.genome.all.supervised.2.Q"
    for k in "${unsupervised_ks[@]}"; do
        "${admixture_exec}" \
            -j"${num_threads}" \
            "onekg.rep_0.genome.all.bed" "${k}"
        cp "onekg.rep_0.genome.all.${k}.Q" \
            "onekg.rep_0.genome.all.unsupervised.${k}.Q"
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
    admixture_q_args+=(
        "--admixture-q-path"
        "${k}=${q_path}"
    )
done

# run empirical-only fastStructure in its Python 2.7 environment.
(
    conda activate "${faststructure_conda_env}"
    for k in "${unsupervised_ks[@]}"; do
        python "${faststructure_structure_py}" \
            -K "${k}" \
            --input="${admixture_prefix}" \
            --output="${faststructure_prefix}" \
            --prior="${faststructure_prior}" \
            --cv="${faststructure_cv}" \
            --seed="${faststructure_seed}"
    done
    python "${faststructure_choose_k_py}" \
        --input="${faststructure_prefix}" \
        | tee "${faststructure_choose_k_path}"
)
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
    faststructure_q_args+=(
        "--faststructure-q-path"
        "${k}=${q_path}"
    )
done
if [[ ! -s "${faststructure_choose_k_path}" ]]; then
    echo "ERROR: missing fastStructure chooseK report" >&2
    exit 1
fi


##### statistics ##############################################################
# parse both ancestry tools against the final FAM order, then combine stats.
python "${project_dir}/job_scripts/python_utils/build_onekg_ancestry.py" \
    --unrels-path "${unrels_path}" \
    --source-fam-path "${source_fam_path}" \
    --admixture-fam-path "${admixture_prefix}.fam" \
    --supervised-q-path "${admixture_prefix}.supervised.2.Q" \
    "${admixture_q_args[@]}" \
    "${faststructure_q_args[@]}" \
    --faststructure-choose-k-path "${faststructure_choose_k_path}" \
    --faststructure-prior "${faststructure_prior}" \
    --faststructure-seed "${faststructure_seed}" \
    --stats-dir "${stats_dir}" \
    --afr-pop "${pops[0]}" \
    --eur-pop "${pops[1]}" \
    --admixed-pop "${pops[2]}"
python "${project_dir}/job_scripts/python_utils/combine_onekg_genome_stats.py" \
    --stats-dir "${stats_dir}" \
    --chroms "${chroms[@]}" \
    --pops "${pops[@]}"
log_msg "done combining 1000 Genomes genome statistics"

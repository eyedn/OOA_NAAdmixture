#!/usr/bin/env bash

###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           calc_sim_genome_rep_stats.sh
###############################################################################

# workflow: merge simulation chromosomes and calculate genome statistics.


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

# load the software environment used for VCF, PLINK, ADMIXTURE, and Python.
module purge
ml gcc/13.3.0 htslib/1.19.1 bcftools/1.19 plink2/2.00a4.3 conda
source /apps/conda/miniforge3/25.3.0/etc/profile.d/conda.sh
conda activate "${conda_env}"

# load shared logging and require one Slurm task per replicate.
project_dir="$(pwd)"
source "${project_dir}/other_scripts/log_msg.sh"

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
# read fixed inputs followed by chromosome and population arrays.
vcf_dir="$1"
plink_bed_dir="$2"
pop_info_dir="$3"
admixture_dir="$4"
king_dir="$5"
stats_dir="$6"
sample_size="$7"
num_reps="$8"
genetic_map="$9"
kin_cutoff="${10}"
admixture_ld_window="${11}"
admixture_ld_step="${12}"
admixture_ld_r2="${13}"
shift 13
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

# derive replicate-specific paths shared by all genome-level stages.
rep="${SLURM_ARRAY_TASK_ID}"
num_threads="${SLURM_CPUS_PER_TASK:-1}"
genome_prefix="${genetic_map}_${rep}_genome_all"
merged_vcf_path="${vcf_dir}/${genome_prefix}.biallelic_snps.vcf.gz"
genome_bed_prefix="${plink_bed_dir}/${genome_prefix}"
sample_metadata_path="${pop_info_dir}/${genetic_map}_${rep}"
sample_metadata_path+="_chr${chroms[0]}.sample_metadata.tsv"
admixture_prefix="${admixture_dir}/${genome_prefix}"
faststructure_prefix="${faststructure_dir}/${genome_prefix}"
faststructure_choose_k_path="${faststructure_prefix}.chooseK.txt"
unrelated_keep_path="${admixture_dir}/${genome_prefix}.king_unrelated.keep"
ld_prune_prefix="${admixture_dir}/${genome_prefix}.ld_prune"
inference_seed="${rep}"

if (( rep < 1 || rep > num_reps )); then
    echo "ERROR: invalid replicate ${rep}; expected 1..${num_reps}" >&2
    exit 1
fi


##### genome PLINK input ######################################################
# create output directories and merge chromosome VCFs for this replicate.
mkdir -p "${admixture_dir}" "${faststructure_dir}" \
    "${king_dir}" "${stats_dir}" "${plink_bed_dir}"

# concatenate chromosome-local VCFs for one autosomal replicate
vcf_paths=()
for chr in "${chroms[@]}"; do
    vcf_paths+=(
        "${vcf_dir}/${genetic_map}_${rep}_chr${chr}_all.biallelic_snps.vcf.gz"
    )
done

log_msg "concatenating chromosome VCFs for rep=${rep}"
bcftools concat \
    --threads "${num_threads}" \
    -Oz \
    -o "${merged_vcf_path}" \
    "${vcf_paths[@]}"
tabix -f -p vcf "${merged_vcf_path}"

log_msg "creating merged autosomal PLINK BED with unique var ids for rep=${rep}"
plink2 \
    --vcf "${merged_vcf_path}" \
    --set-all-var-ids '@:#:$r:$a' \
    --threads "${num_threads}" \
    --make-bed \
    --out "${genome_bed_prefix}"

if [[ ! -s "${genome_bed_prefix}.bed" || ! -s "${genome_bed_prefix}.bim" \
    || ! -s "${genome_bed_prefix}.fam" ]]; then
    echo "ERROR: failed to create merged genome BED set" >&2
    exit 1
fi


##### KING coefficients #######################################################
# run KING for each population and retain unrelated samples for ADMIXTURE.
log_msg "running genome-level KING for rep=${rep}"
: > "${unrelated_keep_path}"
for pop in "${pops[@]}"; do
    subset_path="${king_dir}/${genetic_map}_${rep}_genome_${pop}.subset"
    out_prefix="${king_dir}/${genetic_map}_${rep}_genome_${pop}"
    python "${project_dir}/job_scripts/python_utils/write_sim_pop_subset.py" \
        --subset-path "${subset_path}" \
        --pop "${pop}" \
        --sample-size "${sample_size}" \
        --pops "${pops[@]}"
    plink2 \
        --bfile "${genome_bed_prefix}" \
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

    unrelated_out_prefix="${out_prefix}_unrelated"
    plink2 \
        --bfile "${genome_bed_prefix}" \
        --keep "${retained_path}" \
        --threads "${num_threads}" \
        --make-king \
        --make-king-table \
        --out "${unrelated_out_prefix}"
    cat "${retained_path}" >> "${unrelated_keep_path}"
done

python \
    "${project_dir}/job_scripts/python_utils/write_sim_unrelated_kinship.py" \
    --rep "${rep}" \
    --king-dir "${king_dir}" \
    --stats-dir "${stats_dir}" \
    --genetic-map "${genetic_map}" \
    --pops "${pops[@]}"


##### ADMIXTURE ###############################################################
# LD-prune the combined unrelated samples before supervised ADMIXTURE.
log_msg "LD-pruning genome-level ADMIXTURE samples for rep=${rep}"
plink2 \
    --bfile "${genome_bed_prefix}" \
    --keep "${unrelated_keep_path}" \
    --threads "${num_threads}" \
    --indep-pairwise "${admixture_ld_window}" \
        "${admixture_ld_step}" \
        "${admixture_ld_r2}" \
    --out "${ld_prune_prefix}"

if [[ ! -s "${ld_prune_prefix}.prune.in" ]]; then
    echo "ERROR: missing genome LD-pruned SNP list" >&2
    exit 1
fi

log_msg "creating genome-level ADMIXTURE BED for rep=${rep}"
plink2 \
    --bfile "${genome_bed_prefix}" \
    --keep "${unrelated_keep_path}" \
    --extract "${ld_prune_prefix}.prune.in" \
    --threads "${num_threads}" \
    --make-bed \
    --out "${admixture_prefix}"

log_msg "writing genome-level supervised ADMIXTURE pop file for rep=${rep}"
python "${project_dir}/job_scripts/python_utils/write_sim_admixture_pop.py" \
    --sample-metadata-path "${sample_metadata_path}" \
    --fam-path "${admixture_prefix}.fam" \
    --pop-path "${admixture_prefix}.pop"

# preserve supervised K=2 before unsupervised K=2 reuses ADMIXTURE's default
# output name.
log_msg "running genome-level ADMIXTURE for rep=${rep}"
(
    cd "${admixture_dir}"
    "${admixture_exec}" \
        --supervised \
        -j"${num_threads}" \
        -s "${rep}" \
        "${genome_prefix}.bed" \
        2
    cp "${genome_prefix}.2.Q" "${genome_prefix}.supervised.2.Q"
    for k in "${unsupervised_ks[@]}"; do
        "${admixture_exec}" \
            -j"${num_threads}" \
            -s "${inference_seed}" \
            "${genome_prefix}.bed" \
            "${k}"
        cp "${genome_prefix}.${k}.Q" \
            "${genome_prefix}.unsupervised.${k}.Q"
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

# write fresh genome-wide neutral multi-K tables in final FAM order.
python \
    "${project_dir}/job_scripts/python_utils/build_sim_inference_ancestry.py" \
    --rep "${rep}" \
    --sample-metadata-path "${sample_metadata_path}" \
    --admixture-fam-path "${admixture_prefix}.fam" \
    "${admixture_q_args[@]}" \
    "${faststructure_q_args[@]}" \
    --faststructure-choose-k-path "${faststructure_choose_k_path}" \
    --faststructure-prior "${faststructure_prior}" \
    --faststructure-seed "${inference_seed}" \
    --stats-dir "${stats_dir}"


##### statistics ##############################################################
# aggregate chromosome tables and genome-level KING and ADMIXTURE outputs.
log_msg "aggregating genome-level statistics for rep=${rep}"
python "${project_dir}/job_scripts/python_utils/aggregate_sim_genome_stats.py" \
    --rep "${rep}" \
    --stats-dir "${stats_dir}" \
    --admixture-dir "${admixture_dir}" \
    --king-dir "${king_dir}" \
    --genetic-map "${genetic_map}" \
    --sample-size "${sample_size}" \
    --chroms "${chroms[@]}" \
    --pops "${pops[@]}"

log_msg "done with genome-level statistics for rep=${rep}"

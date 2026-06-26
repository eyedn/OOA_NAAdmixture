#!/usr/bin/env bash

###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology 
#           Mooney Lab
#           ---
#           calc_genome_stats.sh
###############################################################################


#SBATCH --time=1-00:00:00
#SBATCH --partition=qcb
#SBATCH --account=jazlynmo_738
#SBATCH --nodes=1
#SBATCH --output=/home1/karatas/logs/calcGenomeOOANAA/calcGenomeOOANAA.%A_%a.%x.out
#SBATCH --error=/home1/karatas/logs/calcGenomeOOANAA/calcGenomeOOANAA.%A_%a.%x.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=karatas@usc.edu


set -euo pipefail

module purge
ml gcc/13.3.0 htslib/1.19.1 bcftools/1.19 plink2/2.00a4.3 conda
source /apps/conda/miniforge3/25.3.0/etc/profile.d/conda.sh
conda activate OOA_NAAdmixture
export PATH="${HOME}/.conda/envs/OOA_NAAdmixture/bin:${PATH}"

project_dir="$(pwd)"
source "${project_dir}/other_scripts/log_msg.sh"

: "${SLURM_ARRAY_TASK_ID:?ERROR: run as a Slurm array}"


# input variables
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
chroms=()
while [[ "$1" != "--" ]]; do
    chroms+=( "$1" )
    shift
done
shift 1 # skip the "--" from input arguments
pops=( "$@" )

# derived variables
rep="${SLURM_ARRAY_TASK_ID}"
num_threads="${SLURM_CPUS_PER_TASK:-1}"
genome_prefix="${genetic_map}_${rep}_genome_all"
merged_vcf_path="${vcf_dir}/${genome_prefix}.biallelic_snps.vcf.gz"
genome_bed_prefix="${plink_bed_dir}/${genome_prefix}"
sample_metadata_path="${pop_info_dir}/${genetic_map}_${rep}_chr1.sample_metadata.tsv"
admixture_prefix="${admixture_dir}/${genome_prefix}"
unrelated_keep_path="${admixture_dir}/${genome_prefix}.king_unrelated.keep"
ld_prune_prefix="${admixture_dir}/${genome_prefix}.ld_prune"

if (( rep < 1 || rep > num_reps )); then
    echo "ERROR: invalid replicate ${rep}; expected 1..${num_reps}" >&2
    exit 1
fi

mkdir -p "${admixture_dir}" "${king_dir}" "${stats_dir}" "${plink_bed_dir}"

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

log_msg "creating merged autosomal PLINK BED for rep=${rep}"
plink2 \
    --vcf "${merged_vcf_path}" \
    --threads "${num_threads}" \
    --make-bed \
    --out "${genome_bed_prefix}"

if [[ ! -s "${genome_bed_prefix}.bed" || ! -s "${genome_bed_prefix}.bim" \
    || ! -s "${genome_bed_prefix}.fam" ]]; then
    echo "ERROR: failed to create merged genome BED set" >&2
    exit 1
fi

# run genome-level KING, LD pruning, and supervised ADMIXTURE once per replicate
log_msg "running genome-level KING for rep=${rep}"
: > "${unrelated_keep_path}"
for pop in "${pops[@]}"; do
    subset_path="${king_dir}/${genetic_map}_${rep}_genome_${pop}.subset"
    out_prefix="${king_dir}/${genetic_map}_${rep}_genome_${pop}"
    python "${project_dir}/job_scripts/python_utils/write_pop_subset.py" \
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
    cat "${retained_path}" >> "${unrelated_keep_path}"
done

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
python "${project_dir}/job_scripts/python_utils/write_admixture_pop.py" \
    --sample-metadata-path "${sample_metadata_path}" \
    --fam-path "${admixture_prefix}.fam" \
    --pop-path "${admixture_prefix}.pop"

log_msg "running genome-level supervised ADMIXTURE for rep=${rep}"
(
    cd "${admixture_dir}"
    "${HOME}/software/ADMIXTURE/admixture_linux-1.4.0/admixture" \
        --supervised \
        -j"${num_threads}" \
        -s "${rep}" \
        "${genome_prefix}.bed" \
        2
)

log_msg "aggregating genome-level statistics for rep=${rep}"
python "${project_dir}/job_scripts/python_utils/aggregate_genome_stats.py" \
    --rep "${rep}" \
    --stats-dir "${stats_dir}" \
    --admixture-dir "${admixture_dir}" \
    --king-dir "${king_dir}" \
    --genetic-map "${genetic_map}" \
    --chroms "${chroms[@]}" \
    --pops "${pops[@]}"

log_msg "done with genome-level statistics for rep=${rep}"

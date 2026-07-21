#!/usr/bin/env bash

###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           calc_onekg_chr_pop_stats.sh
###############################################################################

# workflow: calculate 1000 Genomes statistics for one chromosome and pop.


##### set up ##################################################################
set -euo pipefail

# load the software environment used for PLINK, VCF, and Python statistics.
module purge
ml gcc/13.3.0 htslib/1.19.1 plink2/2.00a4.3 conda
source /apps/conda/miniforge3/25.3.0/etc/profile.d/conda.sh
conda activate OOA_NAAdmixture

# load shared logging and chromosome-population task mapping.
project_dir="$(pwd)"
source "${project_dir}/other_scripts/log_msg.sh"
source "${project_dir}/other_scripts/map_slurm_task.sh"
: "${SLURM_ARRAY_TASK_ID:?ERROR: run as a Slurm array}"


##### variables ###############################################################
# read fixed inputs followed by chromosome and population arrays.
vcf_prefix="$1"
vcf_suffix="$2"
unrels_path="$3"
fam_path="$4"
chr_lens_path="$5"
intergenic_file="$6"
out_vcf_dir="$7"
out_bed_dir="$8"
king_dir="$9"
stats_dir="${10}"
mutation_rate="${11}"
ld_decay_window_size_bp="${12}"
ld_decay_distance_bin_bp="${13}"
ld_decay_maf_threshold="${14}"
shift 14
shift 1 # skip the "--" from input arguments
chroms=()
while [[ "$1" != "--" ]]; do
    chroms+=( "$1" )
    shift
done
shift 1 # skip the "--" from input arguments
pops=( "$@" )

# map this chromosome-major task ID and derive all output prefixes.
read -r chr pop < <(
    map_slurm_task_to_chr_pop \
        "${SLURM_ARRAY_TASK_ID}" \
        -- "${chroms[@]}" -- "${pops[@]}"
)
num_threads="${SLURM_CPUS_PER_TASK:-1}"
input_vcf="${vcf_prefix}${chr}${vcf_suffix}"
prefix="onekg.rep_0.chr${chr}"
pop_bed_prefix="${out_bed_dir}/${prefix}.${pop}"
common_vcf="${out_vcf_dir}/${prefix}.all.${pop}.vcf.gz"
intergenic_vcf="${out_vcf_dir}/${prefix}.${pop}.intergenic.vcf.gz"
king_prefix="${king_dir}/${prefix}.${pop}"


##### KING coefficients #######################################################
# calculate KING coefficients for unrelated samples in the target population.
mkdir -p "${king_dir}" "${stats_dir}"
log_msg "running unrelated-sample KING chr=${chr} pop=${pop}"
plink2 \
    --bfile "${pop_bed_prefix}" \
    --threads "${num_threads}" \
    --make-king-table \
    --out "${king_prefix}"


##### statistics ##############################################################
# hand shared VCF and KING outputs to the empirical statistics implementation.
python "${project_dir}/job_scripts/python_utils/calc_onekg_stats.py" \
    --analysis-level "chromosome" \
    --vcf-path "${input_vcf}" \
    --ld-vcf-path "${common_vcf}" \
    --intergenic-vcf-path "${intergenic_vcf}" \
    --intergenic-bed-path "${intergenic_file}" \
    --unrels-path "${unrels_path}" \
    --fam-path "${fam_path}" \
    --chr-lens-path "${chr_lens_path}" \
    --king-path "${king_prefix}.kin0" \
    --stats-dir "${stats_dir}" \
    --chrom "${chr}" \
    --pop "${pop}" \
    --mutation-rate "${mutation_rate}" \
    --ld-decay-window-size-bp "${ld_decay_window_size_bp}" \
    --ld-decay-distance-bin-bp "${ld_decay_distance_bin_bp}" \
    --ld-decay-maf-threshold "${ld_decay_maf_threshold}" \
    --pops "${pops[@]}"
log_msg "done with 1000 Genomes statistics chr=${chr} pop=${pop}"

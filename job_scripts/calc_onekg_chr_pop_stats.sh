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
pop_info_dir="$9"
king_dir="${10}"
stats_dir="${11}"
mutation_rate="${12}"
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
common_bed_prefix="${out_bed_dir}/${prefix}.common_maf10.${pop}"
pop_keep="${pop_info_dir}/${prefix}.${pop}.keep"
king_prefix="${king_dir}/${prefix}.${pop}"
ld_prefix="${stats_dir}/${prefix}.${pop}.ld"


##### KING coefficients #######################################################
# calculate KING coefficients for unrelated samples in the target population.
mkdir -p "${king_dir}" "${stats_dir}"
log_msg "running unrelated-sample KING chr=${chr} pop=${pop}"
plink2 \
    --bfile "${pop_bed_prefix}" \
    --threads "${num_threads}" \
    --make-king-table \
    --out "${king_prefix}"


##### linkage disequilibrium ##################################################
# build the common-site PLINK input and calculate unphased LD pairs.
log_msg "creating shared MAF-filtered LD input chr=${chr} pop=${pop}"
plink2 \
    --vcf "${common_vcf}" \
    --set-all-var-ids '@:#:$r:$a' \
    --maf 0.10 \
    --threads "${num_threads}" \
    --make-bed \
    --out "${common_bed_prefix}"

log_msg "ld prune with PLINK chr=${chr} pop=${pop}"
plink2 \
    --bfile "${common_bed_prefix}" \
    --keep "${pop_keep}" \
    --threads "${num_threads}" \
    --indep-pairwise "${admixture_ld_window}" \
        "${admixture_ld_step}" \
        "${admixture_ld_r2}" \
    --out "${ld_prefix}"

#TODO: correct verify step: verify that the LD-pruned SNP list was created.
if [[ ! -s "${ld_prune_prefix}.prune.in" ]]; then
    echo "ERROR: missing LD-pruned SNP list ${ld_prune_prefix}.prune.in" >&2
    exit 1
fi

# accept the PLINK-version-specific LD output suffix that was produced.
ld_path=""
for candidate in "${ld_prefix}.vcor" "${ld_prefix}.vcor1" \
    "${ld_prefix}.ld"; do
    if [[ -s "${candidate}" ]]; then
        ld_path="${candidate}"
        break
    fi
done
if [[ -z "${ld_path}" ]]; then
    echo "ERROR: PLINK did not produce an LD-pair table" >&2
    exit 1
fi


##### statistics ##############################################################
# hand VCF, KING, and LD outputs to the empirical statistics implementation.
python "${project_dir}/job_scripts/python_utils/calc_onekg_stats.py" \
    --analysis-level "chromosome" \
    --vcf-path "${input_vcf}" \
    --intergenic-vcf-path "${intergenic_vcf}" \
    --intergenic-bed-path "${intergenic_file}" \
    --unrels-path "${unrels_path}" \
    --fam-path "${fam_path}" \
    --chr-lens-path "${chr_lens_path}" \
    --king-path "${king_prefix}.kin0" \
    --ld-path "${ld_path}" \
    --stats-dir "${stats_dir}" \
    --chrom "${chr}" \
    --pop "${pop}" \
    --mutation-rate "${mutation_rate}" \
    --pops "${pops[@]}"
log_msg "done with 1000 Genomes statistics chr=${chr} pop=${pop}"

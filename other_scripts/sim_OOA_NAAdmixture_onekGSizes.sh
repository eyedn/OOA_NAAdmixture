#!/usr/bin/env bash

###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology
#           Mooney Lab
#           ---
#           sim_OOA_NAAdmixture_onekGSizes.sh
###############################################################################

#SBATCH --job-name="OOAOnekGSizes"
#SBATCH --array=1-44%100
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=1-00:00:00
#SBATCH --partition=qcb
#SBATCH --account=jazlynmo_738
#SBATCH --nodes=1
#SBATCH --output="/home1/karatas/logs/OOAOnekGSizes/%A_%a.%x.out"
#SBATCH --error="/home1/karatas/logs/OOAOnekGSizes/%A_%a.%x.err"
#SBATCH --mail-type=ALL
#SBATCH --mail-user=karatas@usc.edu

# workflow one rep per autosome for TC and LG, with diploid 1kG sample sizes.


##### set up and task mapping #################################################
set -euo pipefail
: "${SLURM_ARRAY_TASK_ID:?ERROR: run as a Slurm array task}"
if (( SLURM_ARRAY_TASK_ID < 1 || SLURM_ARRAY_TASK_ID > 44 )); then
    echo "ERROR: invalid array task ID: ${SLURM_ARRAY_TASK_ID}" >&2
    exit 1
fi

# Slurm runs a spooled copy; recover the original submitted script path.
script_path="${BASH_SOURCE[0]}"
if [[ "${script_path}" == */slurm_script ]]; then
    script_path=$(scontrol show job "${SLURM_JOB_ID}" -o |
        tr ' ' '\n' | sed -n 's/^Command=//p')
fi
script_dir="$(cd "$(dirname "${script_path}")" && pwd)"
project_dir="$(cd "${script_dir}/.." && pwd)"
source "${script_dir}/log_msg.sh"
chrom=$(((SLURM_ARRAY_TASK_ID - 1) % 22 + 1))
if (( SLURM_ARRAY_TASK_ID <= 22 )); then
    model="TC"
    trajectory_suffix="2T12Consistent"
else
    model="LG"
    trajectory_suffix="largeGrowth"
fi

# select the existing trajectory branch without changing shared configuration.
trajectory_config=$(sed \
    "s|^OUT_DIR=.*|OUT_DIR=/unused/OOA_NAAdmixture_${trajectory_suffix}|" \
    "${script_dir}/const.sh")
# evaluate only the trusted repository configuration with the selected suffix.
eval "${trajectory_config}"
: "${OOA_NAADMIXTURE_CONDA:?ERROR: conda environment is not configured}"
module purge
ml gcc/13.3.0 htslib/1.19.1 bcftools/1.19 conda
source /apps/conda/miniforge3/25.3.0/etc/profile.d/conda.sh
conda activate "${OOA_NAADMIXTURE_CONDA}"
export PATH="${HOME}/.conda/envs/${OOA_NAADMIXTURE_CONDA}/bin:${PATH}"
export PYTHONDONTWRITEBYTECODE=1
export PYTHONPATH="${project_dir}/job_scripts/python_utils:${PYTHONPATH:-}"


##### simulation inputs and output paths ######################################
rep=1
seed=$((1000 * rep + chrom))
out_dir="/home1/karatas/scratch/OOA_NAAdmixture_onekGSizeMatch"
tree_dir="${out_dir}/${model}/trees"
vcf_dir="${out_dir}/${model}/vcfs"
prefix="${GENETIC_MAP}_${rep}_chr${chrom}_all"
mkdir -p "${tree_dir}" "${vcf_dir}"
work_dir=$(mktemp -d "${vcf_dir}/.${prefix}.XXXXXX")
trap 'rm -rf -- "${work_dir}"' EXIT

# pass shared demographic scalars and ordered generation vectors to Python.
export GENERATION_TIME MUTATION_RATE T_AF_YEARS T_OOA_YEARS T_EU0_YEARS
export T_EG_YEARS R_EU0 R_EU R_AF N_A N_AF1 N_B N_EU0 M_AF_B M_AF_EU
export ADMIXTURE_TIME ADMIX_GENERATION_COUNT ADMIX_MIXING_GENERATION_COUNT
export ADMIX_MODERN_GROWTH_RATE CENSUS_TIME_OFFSET GENETIC_MAP MSPRIME_MODEL
admix_ne_csv=$(IFS=,; echo "${ADMIX_NE_BY_GENERATION[*]}")
admix_afr_csv=$(IFS=,; echo "${ADMIX_AFR_PROPS_BY_GENERATION[*]}")
admix_eur_csv=$(IFS=,; echo "${ADMIX_EUR_PROPS_BY_GENERATION[*]}")
admix_prior_csv=$(IFS=,; echo "${ADMIX_PRIORADMIX_PROPS_BY_GENERATION[*]}")


##### simulation and lossless tree verification ################################
log_msg "simulating ${model} chr${chrom}; seed=${seed}; AFR=118 EUR=119 ADX=50"
python - "${chrom}" "${seed}" "${work_dir}/${prefix}" \
    "${admix_ne_csv}" "${admix_afr_csv}" "${admix_eur_csv}" \
    "${admix_prior_csv}" <<'PYTHON'
import os
import sys
import msprime
import stdpopsim
import tszip
from sim_utils.run_simulation import _build_demography

chrom, seed, output_prefix = sys.argv[1:4]
parameters = {
    name: float(os.environ[name.upper()])
    for name in (
        "generation_time", "mutation_rate", "t_af_years", "t_ooa_years",
        "t_eu0_years", "t_eg_years", "r_eu0", "r_eu", "r_af", "n_a",
        "n_af1", "n_b", "n_eu0", "m_af_b", "m_af_eu", "admixture_time",
        "admix_modern_growth_rate", "census_time_offset",
    )
}
for name in ("admix_generation_count", "admix_mixing_generation_count"):
    parameters[name] = int(os.environ[name.upper()])
for name, raw in zip(
    ("admix_ne_by_generation", "admix_afr_props_by_generation",
     "admix_eur_props_by_generation", "admix_prioradmix_props_by_generation"),
    sys.argv[4:8],
):
    parameters[name] = [float(value) for value in raw.split(",")]
demography, _ = _build_demography(**parameters)
population_sizes = {"AFR": 118, "EUR": 119, "ADX": 50}
population_order = [population.name for population in demography.populations]
if population_order != ["AFR", "EUR", "ADX"]:
    raise ValueError(f"Unexpected model population order: {population_order}")
sample_sets = [
    msprime.SampleSet(population_sizes[pop], population=pop, time=0, ploidy=2)
    for pop in population_order
]
sample_names = [
    f"{pop}_{index}"
    for pop in population_order
    for index in range(1, population_sizes[pop] + 1)
]
contig = stdpopsim.get_species("HomSap").get_contig(
    f"chr{chrom}", genetic_map=os.environ["GENETIC_MAP"]
)
ts = msprime.sim_ancestry(
    samples=sample_sets, demography=demography,
    recombination_rate=contig.recombination_map, sequence_length=contig.length,
    ploidy=2, random_seed=int(seed), model=os.environ["MSPRIME_MODEL"],
)
ts = msprime.sim_mutations(
    ts, rate=parameters["mutation_rate"], random_seed=int(seed)
)
if ts.num_samples != 574 or ts.num_individuals != 287:
    raise ValueError("Unexpected tree-sequence sample or individual count")
# check the population of each diploid individual before assigning VCF names.
expected_populations = [
    pop_id for pop_id, pop in enumerate(population_order)
    for _ in range(population_sizes[pop])
]
for individual, pop_id in zip(ts.individuals(), expected_populations):
    if len(individual.nodes) != 2 or any(
        ts.node(node).population != pop_id for node in individual.nodes
    ):
        raise ValueError("Tree-sequence individuals are out of population order")
ts_path = f"{output_prefix}.ts.tsz"
tszip.compress(ts, ts_path)
if not ts.equals(tszip.decompress(ts_path)):
    raise ValueError("Compressed tree sequence failed lossless round-trip")
with open(f"{output_prefix}.vcf", "w", encoding="utf-8") as output:
    ts.write_vcf(
        output, contig_id=f"chr{chrom}", individual_names=sample_names,
        position_transform=lambda positions: [pos + 1 for pos in positions],
    )
PYTHON


##### compress, index, and verify VCF ##########################################
bgzip -f "${work_dir}/${prefix}.vcf"
vcf_path="${work_dir}/${prefix}.vcf.gz"
tabix -f -p vcf "${vcf_path}"
bgzip -t "${vcf_path}"
# read all records as well as the index and compare every ordered sample label.
bcftools view -Ov "${vcf_path}" > /dev/null
tabix -l "${vcf_path}" > /dev/null
bcftools query -l "${vcf_path}" > "${work_dir}/actual_samples"
for pop in AFR EUR ADX; do
    case "${pop}" in
        AFR) count=118 ;;
        EUR) count=119 ;;
        ADX) count=50 ;;
    esac
    for ((index = 1; index <= count; index++)); do
        echo "${pop}_${index}"
    done
done > "${work_dir}/expected_samples"
cmp "${work_dir}/expected_samples" "${work_dir}/actual_samples"
for suffix in ts.tsz vcf.gz vcf.gz.tbi; do
    if [[ ! -s "${work_dir}/${prefix}.${suffix}" ]]; then
        echo "ERROR: missing or empty output: ${prefix}.${suffix}" >&2
        exit 1
    fi
done

# publish only verified compressed trees, VCFs, and indexes.
mv "${work_dir}/${prefix}.ts.tsz" "${tree_dir}/${prefix}.ts.tsz"
mv "${vcf_path}" "${vcf_dir}/${prefix}.vcf.gz"
mv "${vcf_path}.tbi" "${vcf_dir}/${prefix}.vcf.gz.tbi"
log_msg "completed ${model} chr${chrom}; VCF samples=287"

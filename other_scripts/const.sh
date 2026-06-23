#!/usr/bin/env bash

###############################################################################
#           Aydin Loid Karatas
#           ---
#           University of Southern California
#           Department of Quantitative and Computational Biology 
#           Mooney Lab
#           ---
#           const.sh
###############################################################################

# simulation specifications
POPS=(AFR EUR ADX)
NONADMIXED_POP_LABELS=(AFR EUR)
ADMIXED_POP_LABEL="ADX"
GENETIC_MAP="HapMapII_GRCh38"
NUM_REPS="50"
CHR="18"
SAMPLE_SIZE="300"
MSPRIME_MODEL="dtwf"

# output directories
PROJECT_DIR="${HOME}/OOA_NAAdmixture"
OUTDIR="${SCRATCH}/OOA_NAAdmixture_chr${CHR}"
TREE_DIR="${OUTDIR}/trees"
VCF_DIR="${OUTDIR}/vcfs"
PLINK_BED_DIR="${OUTDIR}/plink_beds"
POP_INFO_DIR="${OUTDIR}/pop_info"
ANC_DIR="${OUTDIR}/local_ancestry"
GLOBAL_ANC_DIR="${OUTDIR}/global_ancestry"
ADMIXTURE_DIR="${OUTDIR}/admixture"
KING_DIR="${OUTDIR}/king"
STATS_DIR="${OUTDIR}/stats"

# job array specifications
MAX_JOBS="100"
SIM_CPUS_PER_TASK="4"
SIM_MEM="32G"
STATS_CPUS_PER_TASK="4"
STATS_MEM="64G"

# Tennessen 2012 model specifications adapted from Fu 2013.
GENERATION_TIME=25
MUTATION_RATE="2.36e-8"
T_AF_YEARS=148000
T_OOA_YEARS=51000
T_EU0_YEARS=23000
T_EG_YEARS=5115
R_EU0="0.00307"
R_EU="0.0195"
R_AF="0.0166"
N_A=7310
N_AF1=14474
N_B=1861
N_EU0=1032
M_AF_B="15e-5"
M_AF_EU="2.5e-5"

# admixture specifications inspired from Hacker 2020 and Mooney 2023
ADMIXTURE_TIME=14
CENSUS_TIME_OFFSET="1e-6"
ADMIX_GENERATION_COUNT=15
ADMIX_MIXING_GENERATION_COUNT=10
ADMIX_FOUNDER_AFR_COUNT=183
ADMIX_FOUNDER_EUR_COUNT=25
ADMIX_MODERN_GROWTH_RATE="0.023175"
ADMIX_NE_BY_GENERATION=(
    493.7874
    5755.8703
    15296.0328
    58666.1031
    146967.1421
    260312.8139
    436921.0291
    858109.3961
    1845355.3309
    3322796.2758
    3400701.3175
    3480432.8916
    3562033.8225
    3645547.9384
    3731020.0950
)
ADMIX_AFR_PROPS_BY_GENERATION=(
    0.850000
    0.904820
    0.791384
    0.786692
    0.719992
    0.494960
    0.130456
    0.060000
    0.060000
    0.060000
    0.000000
    0.000000
    0.000000
    0.000000
    0.000000
)
ADMIX_EUR_PROPS_BY_GENERATION=(
    0.150000
    0.080000
    0.080000
    0.080000
    0.080000
    0.080000
    0.080000
    0.030000
    0.030000
    0.030000
    0.000000
    0.000000
    0.000000
    0.000000
    0.000000
)
ADMIX_PRIORADMIX_PROPS_BY_GENERATION=(
    0.000000
    0.015180
    0.128616
    0.133308
    0.200008
    0.425040
    0.789544
    0.910000
    0.910000
    0.910000
    1.000000
    1.000000
    1.000000
    1.000000
    1.000000
)

# statistics specification
KIN_CUTOFF="0.0442"
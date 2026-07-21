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

# shared constants for simulation, empirical data, and Slurm job resources.


##### shared configuration ####################################################
# codebase directory
PROJECT_DIR="/home1/karatas/OOA_NAAdmixture"

# slurm specifications
MAX_JOBS="100"
SIM_CPUS_PER_TASK="4"
SIM_MEM="8G"
STATS_CPUS_PER_TASK="8"
STATS_MEM="64G"
COMB_MEM="64G"
PARTITION="qcbr" # typically use "qcb"; can use "qcbr"
ACCOUNT="qcb_640" # typically use "jazlynmo_738"; can use "qcb_640" 
MAIL_TYPE="ALL"
MAIL_USER="karatas@usc.edu"

# chromosome structure specifications
# CHROMS=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22)
CHROMS=(22)
GENETIC_MAP="HapMapII_GRCh38"

# statistics specification
KIN_CUTOFF="0.0442"
ADMIXTURE_LD_WINDOW="50"
ADMIXTURE_LD_STEP="10"
ADMIXTURE_LD_R2="0.1"
LD_DECAY_WINDOW_SIZE_BP="2000000"
LD_DECAY_DISTANCE_BIN_BP="5000"
LD_DECAY_MAF_THRESHOLD="0.10"


##### simulation configuration ################################################
# simulation labels and settings
POPS=(AFR EUR ADX)
NONADMIXED_POP_LABELS=(AFR EUR)
ADMIXED_POP_LABEL="ADX"
NUM_REPS="50"
SAMPLE_SIZE="500"
MSPRIME_MODEL="dtwf"

# outputs directory for simulations; note OUT_DIR should either ends in a
# "_small" or "_large" to differentiate between the 2 admixed Ne trajectories
OUT_DIR="/home1/karatas/scratch/OOA_NAAdmixture_small"
TREE_DIR="${OUT_DIR}/trees"
PICKLED_DEMO_META="${OUT_DIR}/pickled_demo_meta"
VCF_DIR="${OUT_DIR}/vcfs"
PLINK_BED_DIR="${OUT_DIR}/plink_beds"
POP_INFO_DIR="${OUT_DIR}/pop_info"
ANC_DIR="${OUT_DIR}/local_ancestry"
GLOBAL_ANC_DIR="${OUT_DIR}/global_ancestry"
ADMIXTURE_DIR="${OUT_DIR}/admixture"
KING_DIR="${OUT_DIR}/king"
STATS_DIR="${OUT_DIR}/stats"

# Tennessen 2012 model specifications adapted from Fu 2013 for AFR and EUR
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

# simulation admixture specifications inspired by Hacker 2020 and Mooney 2023
ADMIXTURE_TIME=14
CENSUS_TIME_OFFSET="1e-6"
ADMIX_GENERATION_COUNT=15
ADMIX_MIXING_GENERATION_COUNT=10
ADMIX_FOUNDER_AFR_COUNT=183
ADMIX_FOUNDER_EUR_COUNT=25
ADMIX_MODERN_GROWTH_RATE="0.023175"

if [[ "${OUT_DIR}" == *small ]]; then
    # ADMIX_NE_BY_GENERATION following Tennessen consistency
    ADMIX_NE_BY_GENERATION=(
        493.7874
        5755.8703
        15296.0328
        58666.1031
        146967.1421
        260312.8139
        436921.0291
        447164.9165
        457648.9783
        468378.8455
        479360.2811
        490599.1834
        502101.5888
        513873.6753
        525921.7657
    )
else
    # ADMIX_NE_BY_GENERATION following Hacker/Mooney/Schraiber consistency
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
fi

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


##### empirical configuration #################################################
# empirical labels and settings
ONEKG_POPS=(YRI CEU ASW)
ONEKG_NONADMIXED_POP_LABELS=(YRI CEU)
ONEKG_ADMIXED_POP_LABEL="ASW"

# paths related to GRCH38.p14
INTERGENTIC_FILE="/home1/karatas/references/GRCh38.p14/gencodeV49.canonicalIntergenic.bed"

# input and output directories for the 1000 Genomes empirical analogue
# note, a valid VCF path uses "${ONEKG_VCF_PREFIX}${CHR}${ONEKG_VCF_SUFFIX}"
ONEKG_DIR="/home1/karatas/1000GenomeNYGC_hg38"
ONEKG_UNRELS_FILE="${ONEKG_DIR}/FileInformation/allUnrels_rm3rd_RmSampsGarlicStevenPaper.txt"
ONEKG_VCF_DIR="${ONEKG_DIR}/vcfs_strictMask"
ONEKG_VCF_PREFIX="${ONEKG_VCF_DIR}/CCDG_14151_B01_GRM_WGS_2020-08-05_chr"
ONEKG_VCF_SUFFIX=".filtered.shapeit2-duohmm-phased.nodupmarkers.snps.strict.vcf.gz"
ONEKG_PLINK_BED_DIR="${ONEKG_DIR}/plinkFiles_strictmask"
ONEKG_PLINK_BED_FILE="${ONEKG_PLINK_BED_DIR}/allPops.allChroms.snps.QCIndivsForAuto_UnrelsOnly.bed"
ONEKG_PLINK_BIM_FILE="${ONEKG_PLINK_BED_DIR}/allPops.allChroms.snps.QCIndivsForAuto_UnrelsOnly.bim"
ONEKG_PLINK_FAM_FILE="${ONEKG_PLINK_BED_DIR}/allPops.allChroms.snps.QCIndivsForAuto_UnrelsOnly.fam"
ONEKG_CHR_LENS="/home1/karatas/proj/1000GenomeNYGC_hg38_karatas/ONEKG_chr_lens.tsv"

ONEKG_OUT_DIR="/home1/karatas/scratch/OOA_NAAdmixture_1kG"
ONEKG_OUT_VCF_DIR="${ONEKG_OUT_DIR}/vcfs"
ONEKG_OUT_PLINK_BED_DIR="${ONEKG_OUT_DIR}/plink_beds"
ONEKG_OUT_POP_INFO_DIR="${ONEKG_OUT_DIR}/pop_info"
# ONEKG_OUT_ANC_DIR="${ONEKG_OUT_DIR}/local_ancestry"
# ONEKG_OUT_GLOBAL_ANC_DIR="${ONEKG_OUT_DIR}/global_ancestry"
ONEKG_OUT_ADMIXTURE_DIR="${ONEKG_OUT_DIR}/admixture"
ONEKG_OUT_KING_DIR="${ONEKG_OUT_DIR}/king"
ONEKG_OUT_STATS_DIR="${ONEKG_OUT_DIR}/stats"

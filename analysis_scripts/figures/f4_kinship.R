# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology
# Mooney Lab
# ___
# f4_kinship.R
# ______________________________________________________________________________


# set up ----
setwd("~/OOA_NAAdmixture/analysis_scripts/")
library(tidyverse)
library(patchwork)
if (
  !exists("kinship.bootstrap.tcd.1kg.datatype.interval", envir = .GlobalEnv)
) {
  source("kinship.R")
}


# main ----
f4.kinship <- kinship.bootstrap.tcd.1kg.datatype.interval +
  labs(title = NULL, subtitle = NULL) +
  theme(strip.text = element_blank())
print(f4.kinship)


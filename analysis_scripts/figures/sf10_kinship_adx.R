# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology
# Mooney Lab
# ___
# sf10_kinship_adx.R
# ______________________________________________________________________________


# set up ----
setwd("~/OOA_NAAdmixture/analysis_scripts/")
library(tidyverse)
library(patchwork)
if (
  !exists("kinship.bootstrap.all.datatypes.adx.asw", envir = .GlobalEnv)
) {
  source("kinship.R")
}


# main ----
sf10.kinship.adx <- kinship.bootstrap.all.datatypes.adx.asw +
  labs(title = NULL, subtitle = NULL) +
  theme(strip.text = element_blank())
print(sf10.kinship.adx)


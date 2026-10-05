# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology
# Mooney Lab
# ___
# st2_mixing_proportions.R
# ______________________________________________________________________________


# set up ----
setwd("~/OOA_NAAdmixture/analysis_scripts/")
library(tidyverse)
OUTPUT.DIR <- "/home1/karatas/proj/OOA_NAAdmixture_data"
if (
  !exists("AA.mixing.props", envir = .GlobalEnv)
) {
  source("calc_ADX_Ne.R")
}


# main ----
st2.mixing.proportions <- AA.mixing.props %>%
  rename(generation = g)
readr::write_csv(
  st2.mixing.proportions,
  file.path(OUTPUT.DIR, "st2_mixing_proportions.csv")
)
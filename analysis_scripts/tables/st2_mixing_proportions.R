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
  !exists("admix.tbl", envir = .GlobalEnv)
) {
  source("calc_ADX_Ne.R")
}


# main ----
st2.mixing.proportions <- admix.tbl %>%
  transmute(
    generation,
    s.afr = afr,
    s.eur = eur,
    h = prior.admix
  )
readr::write_csv(
  st2.mixing.proportions,
  file.path(OUTPUT.DIR, "st2_mixing_proportions.csv")
)

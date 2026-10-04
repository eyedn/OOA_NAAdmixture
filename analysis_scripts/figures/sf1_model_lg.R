# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology
# Mooney Lab
# ___
# sf1_model_lg.R
# ______________________________________________________________________________


# set up ----
setwd("~/OOA_NAAdmixture/analysis_scripts/")
library(tidyverse)
library(patchwork)
OUTPUT.DIR <- "/home1/karatas/proj/OOA_NAAdmixture_data"
if (!exists("lg.combined.ne.plot", envir = .GlobalEnv)) {
  source("calc_ADX_Ne.R")
}


# main ----
sf1.model.lg <- lg.combined.ne.plot +
  labs(title = NULL, subtitle = NULL) +
  guides(linetype = "none")

dir.create(OUTPUT.DIR, recursive = TRUE, showWarnings = FALSE)
saveRDS(sf1.model.lg, file.path(OUTPUT.DIR, "sf1_model_lg.rds"))
print(sf1.model.lg)

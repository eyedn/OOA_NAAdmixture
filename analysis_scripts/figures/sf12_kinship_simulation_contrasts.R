# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology
# Mooney Lab
# ___
# sf12_kinship_simulation_contrasts.R
# ______________________________________________________________________________


# set up ----
setwd("~/OOA_NAAdmixture/analysis_scripts/")
library(tidyverse)
library(patchwork)
OUTPUT.DIR <- "/home1/karatas/proj/OOA_NAAdmixture_data"
if (
  !exists("kinship.simulation.contrast.plots", envir = .GlobalEnv)
) {
  source("kinship.R")
}


# main ----
sf12.kinship.simulation.contrasts <- 
  kinship.simulation.contrast.plots$bonferroni +
  labs(title = NULL, subtitle = NULL) +
  theme(legend.position = "none")
dir.create(OUTPUT.DIR, recursive = TRUE, showWarnings = FALSE)
saveRDS(
  sf12.kinship.simulation.contrasts,
  file.path(OUTPUT.DIR, "sf12_kinship_simulation_contrasts.rds")
)
print(sf12.kinship.simulation.contrasts)

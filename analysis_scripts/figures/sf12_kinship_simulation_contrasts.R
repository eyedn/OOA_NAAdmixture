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
print(sf12.kinship.simulation.contrasts)


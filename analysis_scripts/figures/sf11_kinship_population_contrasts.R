# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology
# Mooney Lab
# ___
# sf11_kinship_population_contrasts.R
# ______________________________________________________________________________


# set up ----
setwd("~/OOA_NAAdmixture/analysis_scripts/")
library(tidyverse)
library(patchwork)
OUTPUT.DIR <- "/home1/karatas/proj/OOA_NAAdmixture_data"
if (
  !exists("kinship.population.contrast.plots", envir = .GlobalEnv) ||
  !exists("kinship.empirical.contrast.plots", envir = .GlobalEnv)
) {
  source("kinship.R")
}


# main ----
sf11.kinship.population.contrasts <- 
  (
    kinship.population.contrast.plots$bonferroni +
      labs(title = NULL, subtitle = NULL) +
      theme(axis.title.x = element_blank(), legend.position = "none")
    ) / (
      kinship.empirical.contrast.plots$bonferroni +
        labs(title = NULL, subtitle = NULL) +
        theme(legend.position = "none")
      ) +
  plot_annotation(tag_levels = 'A', tag_suffix = '.)') +
  plot_layout(heights = c(2, 1))
dir.create(OUTPUT.DIR, recursive = TRUE, showWarnings = FALSE)
saveRDS(
  sf11.kinship.population.contrasts,
  file.path(OUTPUT.DIR, "sf11_kinship_population_contrasts.rds")
)
print(sf11.kinship.population.contrasts)

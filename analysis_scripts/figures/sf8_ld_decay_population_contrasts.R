# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology
# Mooney Lab
# ___
# sf8_ld.decay_population_contrasts.R
# ______________________________________________________________________________


# set up ----
setwd("~/OOA_NAAdmixture/analysis_scripts/")
library(tidyverse)
library(patchwork)
library(ggh4x)
OUTPUT.DIR <- "/home1/karatas/proj/OOA_NAAdmixture_data"
if (
  !exists("ld.population.contrast.plots", envir = .GlobalEnv) || 
  !exists("ld.population.contrast.plots", envir = .GlobalEnv)
) {
  source("diversity.R")
}


# main ----
sf8.ld.decay.population.contrasts <- 
  (
    ld.population.contrast.plots$tc.tcd$bonferroni +
      labs(title = NULL, subtitle = NULL) +
      theme(axis.title.x = element_blank())
    ) / (
      ld.population.contrast.plots$empirical +
        labs(title = NULL, subtitle = NULL)
      ) +
  plot_annotation(tag_levels = 'A', tag_suffix = '.)') +
  plot_layout(heights = c(2, 1))
  
dir.create(OUTPUT.DIR, recursive = TRUE, showWarnings = FALSE)
saveRDS(
  sf8.ld.decay.population.contrasts,
  file.path(OUTPUT.DIR, "sf8_ld_decay_population_contrasts.rds")
)
print(sf8.ld.decay.population.contrasts)

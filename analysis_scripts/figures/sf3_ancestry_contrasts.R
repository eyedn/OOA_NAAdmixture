# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology
# Mooney Lab
# ___
# sf3_ancestry_contrasts.R
# ______________________________________________________________________________


# set up ----
setwd("~/OOA_NAAdmixture/analysis_scripts/")
library(tidyverse)
library(patchwork)
OUTPUT.DIR <- "/home1/karatas/proj/OOA_NAAdmixture_data"
if (
  !exists(
    "ancestry.bootstrap.empirical.bonferroni.genome.comparisons", 
    envir = .GlobalEnv
    ) || 
  !exists(
    "ancestry.bootstrap.simulation.bonferroni.comparisons", envir = .GlobalEnv
  )
) {
  source("ancestry.R")
}


# main ----
sf3.ancestry.contrasts <- 
  (ancestry.bootstrap.empirical.bonferroni.genome.comparisons +
     labs(title = NULL, subtitle = NULL) +
     facet_grid(
       cols = vars(contrast), rows = vars(statistic), scales = "free_y", 
       labeller = labeller(
         contrast = ANCESTRY.PLOT.STYLES$contrast.labels,
         statistic = c(mean = "Mean", sd = "SD"))
     ) +
     theme(legend.position = "none", axis.title.x = element_blank())
  ) / (
    ancestry.bootstrap.simulation.bonferroni.comparisons +
      labs(title = NULL, subtitle = NULL) +
      facet_grid(
        cols = vars(contrast), rows = vars(statistic), scales = "free_y", 
        labeller = labeller(
          contrast = ANCESTRY.PLOT.STYLES$contrast.labels,
          statistic = c(mean = "Mean", sd = "SD"))
      ) +
      theme(legend.position = "none")
    ) +
  plot_annotation(tag_levels = 'A', tag_suffix = '.)')
dir.create(OUTPUT.DIR, recursive = TRUE, showWarnings = FALSE)
saveRDS(
  sf3.ancestry.contrasts,
  file.path(OUTPUT.DIR, "sf3_ancestry_contrasts.rds")
)
print(sf3.ancestry.contrasts)

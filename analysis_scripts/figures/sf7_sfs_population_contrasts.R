# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology
# Mooney Lab
# ___
# sf7_sfs_population_contrasts.R
# ______________________________________________________________________________


# set up ----
setwd("~/OOA_NAAdmixture/analysis_scripts/")
library(tidyverse)
library(patchwork)
library(ggh4x)
OUTPUT.DIR <- "/home1/karatas/proj/OOA_NAAdmixture_data"
if (
  !exists("sfs.population.contrast.plots", envir = .GlobalEnv) || 
  !exists("sfs.population.contrast.plots", envir = .GlobalEnv)
) {
  source("sfs")
}


# main ----
sf7.sfs.population.contrasts <- 
  (sfs.population.contrast.plots$tc.tcd$bonferroni + 
     labs(title = NULL, subtitle = NULL) +
     facet_grid(
       interaction(data.type, measure, sep = " ") ~ contrast, 
       scales = "free_y",
       labeller = labeller(
         contrast = PLOT.STYLES$contrast.labels,
         .rows = as_labeller(c(
           "Simulation_2T12Consistent count" = "T.C. count",
           "Simulation_2T12Consistent_simDown count" = "T.C.D. count",
           "Simulation_2T12Consistent proportion" = "T.C. prop.",
           "Simulation_2T12Consistent_simDown proportion" = "T.C.D. prop."
         ))
       )
     ) +
     theme(
       legend.position = "none", axis.title.x = element_blank())
  ) / (
    sfs.population.contrast.plots$empirical + 
      labs(title = NULL, subtitle = NULL) +
      facet_grid(
        cols = vars(contrast), rows = vars(measure), scales = "free_y", 
        labeller = labeller(
          contrast = SFS.CONTRAST.LABELS,
          measure = c("count" = "count", "proportion" = "prop.")
        )
      ) +
      theme(legend.position = "none")
  ) +
  plot_annotation(tag_levels = 'A', tag_suffix = '.)') +
  plot_layout(heights = c(2, 1))
dir.create(OUTPUT.DIR, recursive = TRUE, showWarnings = FALSE)
saveRDS(
  sf7.sfs.population.contrasts,
  file.path(OUTPUT.DIR, "sf7_sfs_population_contrasts.rds")
)
print(sf7.sfs.population.contrasts)

# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology
# Mooney Lab
# ___
# sf6_diversity_population_contrasts.R
# ______________________________________________________________________________


# set up ----
setwd("~/OOA_NAAdmixture/analysis_scripts/")
library(tidyverse)
library(patchwork)
library(ggh4x)
OUTPUT.DIR <- "/home1/karatas/proj/OOA_NAAdmixture_data"
if (
  !exists("diversity.population.contrast.plots", envir = .GlobalEnv) || 
  !exists("diversity.population.contrast.plots", envir = .GlobalEnv)
) {
  source("diversity.R")
}


# main ----
sf6.diversity.population.contrasts <- 
  (diversity.population.contrast.plots$tc.tcd$bonferroni +
     labs(title = NULL, subtitle = NULL) +
     facet_grid(
       interaction(data.type, stat, sep = " ") ~ contrast, 
       scales = "free_y",
       labeller = labeller(
         contrast = PLOT.STYLES$contrast.labels,
         .rows = as_labeller(c(
           "Simulation_2T12Consistent pi" = "T.C. π",
           "Simulation_2T12Consistent_simDown pi" = "T.C.D. π",
           "Simulation_2T12Consistent theta" = "T.C. θ",
           "Simulation_2T12Consistent_simDown theta" = "T.C.D. θ"
           ))
         )
     ) +
     theme(
       legend.position = "none", axis.title.x = element_blank(),
       )
  ) / (
    diversity.population.contrast.plots$empirical +
      labs(title = NULL, subtitle = NULL) +
      facet_grid(
        cols = vars(contrast), rows = vars(stat), scales = "free_y", 
        labeller = labeller(
          contrast = PLOT.STYLES$contrast.labels.emp,
          stat = c(pi = "π", theta = "θ"))
      ) +
      theme(
        legend.position = "none"
        ) 
  ) +
  plot_annotation(tag_levels = 'A', tag_suffix = '.)') +
  plot_layout(heights = c(2.5, 1))
dir.create(OUTPUT.DIR, recursive = TRUE, showWarnings = FALSE)
saveRDS(
  sf6.diversity.population.contrasts,
  file.path(OUTPUT.DIR, "sf6_diversity_population_contrasts.rds")
)
print(sf6.diversity.population.contrasts)

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


sf3.contrast.labels <- c(
  `TC-Emp` = "T.C. ADX - ASW", `TCD-Emp` = "T.C.D. ADX- ASW",
  `LG-Emp` = "L.G. ADX - ASW", `LGD-Emp` = "L.G.D. ADX - ASW",
  `TC-TCD` = "T.C. ADX - T.C.D. ADX", `LG-LGD` = "L.G. ADX - L.G.D. ADX",
  `TC-LG` = "T.C. ADX - L.G. ADX"
  )

# main ----
sf3.ancestry.contrasts <- 
  (ancestry.bootstrap.empirical.bonferroni.genome.comparisons +
     labs(title = NULL, subtitle = NULL) +
     facet_grid(
       cols = vars(contrast), rows = vars(statistic), scales = "free_y", 
       labeller = labeller(
         contrast = sf3.contrast.labels,
         statistic = c(mean = "Mean", sd = "SD"))
     ) +
     theme(legend.position = "none", axis.title.x = element_blank())
  ) / (
    ancestry.bootstrap.simulation.bonferroni.comparisons +
      labs(title = NULL, subtitle = NULL) +
      facet_grid(
        cols = vars(contrast), rows = vars(statistic), scales = "free_y", 
        labeller = labeller(
          contrast = sf3.contrast.labels,
          statistic = c(mean = "Mean", sd = "SD"))
      ) +
      theme(legend.position = "none")
    ) +
  plot_annotation(tag_levels = 'A', tag_suffix = '.)')
print(sf3.ancestry.contrasts)

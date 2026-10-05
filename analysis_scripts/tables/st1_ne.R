# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology
# Mooney Lab
# ___
# st1_ne.R
# ______________________________________________________________________________


# set up ----
setwd("~/OOA_NAAdmixture/analysis_scripts/")
library(tidyverse)
OUTPUT.DIR <- "/home1/karatas/proj/OOA_NAAdmixture_data"
if (
  !exists("AA.ne", envir = .GlobalEnv) ||
  !exists("LG.ne", envir = .GlobalEnv)
    ) {
  source("calc_ADX_Ne.R")
}


# main ----
TC <- AA.ne %>%
  select(
    generation, generation.start.year, generation.end.year, pop.end, imported, 
    births, birth.rate, admix.ne, import.ne, ne
    ) %>%
  rename(admix.ne.tc = admix.ne, import.ne.tc = import.ne, ne.tc = ne)

LG <- LG.ne %>%
  select(
    generation, generation.start.year, generation.end.year, pop.end, imported, 
    births, birth.rate, admix.ne, import.ne, ne
    ) %>%
  rename(admix.ne.lg = admix.ne, import.ne.lg = import.ne, ne.lg = ne)

st1.ne <- TC %>% left_join(
  LG, by = c(
    "generation", "generation.start.year", "generation.end.year", "pop.end",
    "imported", "births", "birth.rate"
    )
  )
readr::write_csv(
  st1.ne,
  file.path(OUTPUT.DIR, "st1_ne.csv")
)

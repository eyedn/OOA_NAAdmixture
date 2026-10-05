# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology
# Mooney Lab
# ___
# st5_diversity_pi.R
# ______________________________________________________________________________


# set up ----
setwd("~/OOA_NAAdmixture/analysis_scripts/")
library(tidyverse)
OUTPUT.DIR <- "/home1/karatas/proj/OOA_NAAdmixture_data"
if (
  !exists("diversity.simulation.summary", envir = .GlobalEnv) ||
  !exists("diversity.emp.intergenic.chromosome", envir = .GlobalEnv) ||
  !exists("diversity.emp.intergenic.genome", envir = .GlobalEnv)
) {
  source("diversity.R")
}
chrom.levels <- c("all", as.character(1:22))


# main ----
pi.sim.chrom.data <- diversity.simulation.summary %>%
  filter((mask == "Intergenic" | is.na(mask)), stat == "pi") %>%
  select(data.type, pop, chrom, mean) %>%
  pivot_wider(names_from = chrom, values_from = mean)

pi.emp.chrom.data <- diversity.emp.intergenic.chromosome %>%
  filter(mask == "Intergenic", stat == "pi") %>%
  select(data.type, pop, chrom, value) %>%
  pivot_wider(names_from = chrom, values_from = value)

pi.emp.genome.data <- diversity.emp.intergenic.genome %>%
  filter(mask == "Intergenic", stat == "pi") %>%
  select(data.type, pop, chrom, value) %>%
  pivot_wider(names_from = chrom, values_from = value)

pi.emp.data <- full_join(pi.emp.chrom.data, pi.emp.genome.data)

st5.diversity.pi <- bind_rows(
  pi.sim.chrom.data, 
  pi.emp.data
  ) %>%
  select(
    any_of(names(.)[!names(.) %in% chrom.levels]), any_of(chrom.levels)
    ) %>%
  rename(`genome-wide` = all)
readr::write_csv(
  st5.diversity.pi,
  file.path(OUTPUT.DIR, "st5_diversity_pi.csv")
)
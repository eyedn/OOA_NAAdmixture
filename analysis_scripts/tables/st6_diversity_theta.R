# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology
# Mooney Lab
# ___
# st6_diversity_theta.R
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
theta.sim.chrom.data <- diversity.simulation.summary %>%
  filter((mask == "Intergenic" | is.na(mask)), stat == "theta") %>%
  select(data.type, pop, chrom, mean) %>%
  pivot_wider(names_from = chrom, values_from = mean)

theta.emp.chrom.data <- diversity.emp.intergenic.chromosome %>%
  filter(mask == "Intergenic", stat == "theta") %>%
  select(data.type, pop, chrom, value) %>%
  pivot_wider(names_from = chrom, values_from = value)

theta.emp.genome.data <- diversity.emp.intergenic.genome %>%
  filter(mask == "Intergenic", stat == "theta") %>%
  select(data.type, pop, chrom, value) %>%
  pivot_wider(names_from = chrom, values_from = value)

theta.emp.data <- full_join(theta.emp.chrom.data, theta.emp.genome.data)

st6.diversity.theta <- bind_rows(
  theta.sim.chrom.data, 
  theta.emp.data
) %>%
  select(
    any_of(names(.)[!names(.) %in% chrom.levels]), any_of(chrom.levels)
  ) %>%
  rename(`genome-wide` = all)
readr::write_csv(
  st6.diversity.theta,
  file.path(OUTPUT.DIR, "st6_diversity_theta.csv")
)
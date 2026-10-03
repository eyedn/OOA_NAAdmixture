# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology
# Mooney Lab
# ___
# generate_all.R
# ______________________________________________________________________________


# set up ----
setwd("~/OOA_NAAdmixture/analysis_scripts/")
figure.scripts <- list.files(
  "figures",
  pattern = "\\.[Rr]$",
  full.names = TRUE
)


# main ----
lapply(sort(figure.scripts), source)

# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology
# Mooney Lab
# ___
# tmp_kinship_histogram.R
# ______________________________________________________________________________

# pattern: Imperative Shell


# reuse setup and helpers without running the main analysis ----
local({
  setwd("~/OOA_NAAdmixture/analysis_scripts/")
  kinship.lines <- readLines("kinship.R", warn = FALSE)
  analysis.start <- match("# analysis ----", kinship.lines)
  if (is.na(analysis.start)) {
    stop("Could not locate the analysis section in kinship.R")
    }
  eval(parse(text = kinship.lines[seq_len(analysis.start - 1L)]),
    envir = environment())

  KINSHIP.BIN.WIDTH <- 0.1
  TEMP.OUTPUT.DIR <- file.path(KINSHIP.OUTPUT.DIR, "tmp_kinship_histogram")
  TEMP.OUTPUT.STEM <- paste0(
    "kinship.bootstrap.tcd.1kg.datatype.interval.bin0.1.ncol1"
    )


  # rebuild the original common-bin summary with the temporary width ----
  kinship.directories <- c(
    KINSHIP.SIM.TC.DATA.DIR, KINSHIP.SIMDOWN.TC.DATA.DIR,
    KINSHIP.SIM.LG.DATA.DIR, KINSHIP.SIMDOWN.LG.DATA.DIR,
    KINSHIP.EMPIRICAL.DATA.DIR
    )
  kinship.data <- map2_dfr(
    kinship.directories, KINSHIP.SOURCE.LEVELS,
    function(directory, source) {
      read.kinship.chromosomes(
        directory, KINSHIP.SELECTED.CHROMOSOMES, source
        )
      }
    )
  kinship.breaks <- make.kinship.breaks(kinship.data, KINSHIP.BIN.WIDTH)
  kinship.summary <- summarize.bootstrap.kinship(kinship.data, kinship.breaks)


  # retain the original display window and stack source facets vertically ----
  kinship.bootstrap.tcd.1kg.datatype.interval <- make.bootstrap.kinship.plot(
    kinship.summary, kinship.breaks,
    KINSHIP.PLOT.STYLES$plot.configs$tcd.1kg,
    "Pairwise KING kinship: chromosome 1 TCD and 1kG",
    x.limits = c(-0.1, 0.05), tag = "tcd.1kg",
    facet.by.datatype = FALSE
    ) +
    facet_wrap(
      ~data.type, ncol = 1,
      labeller = labeller(data.type = KINSHIP.SOURCE.LABELS)
      )


  # save only the temporary figure ----
  dir.create(TEMP.OUTPUT.DIR, recursive = TRUE, showWarnings = FALSE)
  saveRDS(
    kinship.bootstrap.tcd.1kg.datatype.interval,
    file.path(TEMP.OUTPUT.DIR, paste0(TEMP.OUTPUT.STEM, ".rds"))
    )
  ggsave(
    file.path(TEMP.OUTPUT.DIR, paste0(TEMP.OUTPUT.STEM, ".pdf")),
    plot = kinship.bootstrap.tcd.1kg.datatype.interval,
    width = 10, height = 10, units = "in"
    )
  print(kinship.bootstrap.tcd.1kg.datatype.interval)
  message("Temporary histogram saved in: ", TEMP.OUTPUT.DIR)
  })

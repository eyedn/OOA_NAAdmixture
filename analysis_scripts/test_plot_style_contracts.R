# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology
# Mooney Lab
# ___
# test_plot_style_contracts.R
# ______________________________________________________________________________


options(scipen = 999)

# pattern: Mixed (unavoidable)
# Reason: This synthetic harness loads analysis helpers and renders plot objects.

library(tidyverse)


# load a script's constants and helpers without its production analysis section.
load.plot.helpers <- function(path) {
  lines <- readLines(path, warn = FALSE)
  analysis.start <- grep("^# analysis ----", lines)[1L]
  lines <- lines[seq_len(analysis.start - 1L)]
  lines <- lines[!grepl("^library\\(", lines)]
  environment <- new.env(parent = globalenv())
  eval(parse(text = lines), envir = environment)
  return(environment)
  }


# assert a single style-contract condition with an informative message.
assert.style <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
  }


# return a simple synthetic comparison data frame used by discrete contrast plots.
make.discrete.contrast.data <- function(contrasts, x.column) {
  data <- tibble(
    contrast = factor(contrasts, levels = contrasts),
    difference = c(0.20, -0.10, 0.05)[seq_along(contrasts)],
    ci.lower = c(0.05, -0.20, -0.10)[seq_along(contrasts)],
    ci.upper = c(0.30, 0.01, 0.20)[seq_along(contrasts)],
    point.color = c("red", as.character(contrasts[-1L]))
    )
  data[[x.column]] <- 1L
  return(data)
  }


# validate all scripts provide the shared pure plot-theme helper.
script.paths <- c(
  "analysis_scripts/ancestry.R",
  "analysis_scripts/diversity.R",
  "analysis_scripts/kinship.R",
  "analysis_scripts/ld_decay.R",
  "analysis_scripts/sfs.R"
  )
assert.style(
  all(vapply(script.paths, function(path) {
    any(grepl("^options\\(scipen = 999\\)$", readLines(path)))
    }, logical(1))),
  "Not all production scripts disable scientific notation"
  )
helpers <- lapply(script.paths, load.plot.helpers)
names(helpers) <- basename(script.paths)
for (name in names(helpers)) {
  environment <- helpers[[name]]
  assert.style(
    exists("apply.standard.plot.theme", envir = environment),
    paste(name, "does not define apply.standard.plot.theme")
    )
  plot <- ggplot(tibble(x = 1, y = 1), aes(x, y)) + geom_point()
  themed <- environment$apply.standard.plot.theme(plot)
  assert.style(
    identical(themed$theme$panel.grid.minor, element_blank()),
    paste(name, "does not remove minor gridlines")
    )
  assert.style(
    identical(themed$theme$strip.background, element_blank()),
    paste(name, "does not use transparent facet strips")
    )
  assert.style(
    identical(themed$theme$strip.text$face, "plain"),
    paste(name, "does not use plain facet labels")
    )
  }


# verify palette contracts and method-aware ancestry metadata.
cyan.palette <- c("#BFFBFF", "#16ACBD", "#00606F")
assert.style(
  identical(unname(helpers$diversity.R$PLOT.STYLES$contrast.colors), cyan.palette),
  "Diversity population contrasts do not use the cyan palette"
  )
assert.style(
  identical(unname(helpers$kinship.R$KINSHIP.CONTRAST.COLORS[1:3]), cyan.palette),
  "Kinship population contrasts do not use the cyan palette"
  )
assert.style(
  identical(unname(helpers$sfs.R$SFS.CONTRAST.COLORS[1:3]), cyan.palette),
  "SFS population contrasts do not use the cyan palette"
  )
assert.style(
  identical(
    unname(helpers$ancestry.R$PLOT.STYLES$empirical.colors),
    rep("#E44B8D", 2L)
    ),
  "Ancestry empirical inference colors are not consistently ASW pink"
  )
assert.style(
  identical(
    helpers$ancestry.R$ancestry.output.filename("plot", "rds"),
    "plot.admixture.rds"
    ),
  "Ancestry output filenames are not method tagged"
  )
assert.style(
  identical(
    helpers$ancestry.R$ancestry.inference.subtitle(),
    "Empirical inference: ADMIXTURE"
    ),
  "Ancestry subtitle does not identify the inference method"
  )
helpers$ancestry.R$PLOT.EMPIRICAL.METHOD <- "fastStructure"
assert.style(
  identical(
    helpers$ancestry.R$ancestry.output.filename("plot", "rds"),
    "plot.faststructure.rds"
    ),
  "fastStructure ancestry outputs are not method tagged"
  )
helpers$ancestry.R$PLOT.EMPIRICAL.METHOD <- "ADMIXTURE"


# return the first built layer that uses a requested geom class.
built.layer <- function(plot, geom.class) {
  index <- which(vapply(
    plot$layers,
    function(layer) inherits(layer$geom, geom.class), logical(1)
    ))[1L]
  return(ggplot_build(plot)$data[[index]])
  }


# validate shared dot, fill, outline, and interval styling.
assert.discrete.dot.style <- function(plot, contrast.colors, label) {
  point.layers <- vapply(plot$layers, function(layer) {
    inherits(layer$geom, "GeomPoint") &&
      identical(layer$aes_params$shape, 21)
    }, logical(1))
  assert.style(any(point.layers), paste(label, "do not use shape-21 dots"))
  points <- built.layer(plot, "GeomPoint")
  assert.style(
    all(unname(contrast.colors) %in% unique(points$fill)),
    paste(label, "do not use cyan contrast fills")
    )
  assert.style(
    all(c("red", unname(contrast.colors)[2L]) %in% unique(points$colour)),
    paste(label, "do not encode significant and non-significant outlines")
    )
  assert.style(
    all(built.layer(plot, "GeomErrorbar")$colour == "black"),
    paste(label, "do not use black confidence intervals")
    )
  assert.style(
    all(built.layer(plot, "GeomErrorbar")$width == 0),
    paste(label, "do not use zero-width confidence caps")
    )
  }


# verify discrete contrast plots use filled dots and black interval CIs.
ancestry <- helpers$ancestry.R
ancestry.contrasts <- c("TC-TCD", "LG-LGD", "TC-LG")
ancestry.data <- make.discrete.contrast.data(ancestry.contrasts, "plot.x") %>%
  mutate(
    chromosome = factor("1", levels = "1"),
    statistic = factor("mean", levels = "mean"),
    facet.label = factor(NA_character_)
    )
ancestry.plot <- ancestry$make.ancestry.bootstrap.contrast.plot(
  ancestry.data, "Synthetic ancestry", facet.statistics = FALSE
  )
assert.discrete.dot.style(
  ancestry.plot, ancestry$PLOT.STYLES$contrast.colors[ancestry.contrasts],
  "Ancestry contrasts"
  )

diversity <- helpers$diversity.R
diversity.contrasts <- c("AFR-ADX", "AFR-EUR", "ADX-EUR")
diversity.data <- make.discrete.contrast.data(diversity.contrasts, "plot.x") %>%
  mutate(
    chrom = factor("1", levels = "1"),
    data.type = factor("Simulation_2T12Consistent",
      levels = diversity$SOURCE.LEVELS),
    stat = factor("pi", levels = c("pi", "theta"))
    )
diversity.plot <- diversity$make.diversity.contrast.plot(
  list(simulation = diversity.data), "95", "population"
  )
assert.discrete.dot.style(
  diversity.plot, diversity$PLOT.STYLES$contrast.colors,
  "Diversity contrasts"
  )
assert.style(
  all(built.layer(ancestry.plot, "GeomHline")$linetype == "dashed") &&
    all(built.layer(ancestry.plot, "GeomHline")$linewidth == 1),
  "Ancestry contrast zero baseline does not use the standard dashed stroke"
  )

kinship <- helpers$kinship.R
kinship.contrasts <- c("AFR-ADX", "AFR-EUR", "ADX-EUR")
kinship.data <- make.discrete.contrast.data(kinship.contrasts, "xmid") %>%
  mutate(xmid = c(-0.10, -0.09, -0.08))
kinship.plot <- kinship$make.kinship.contrast.plot(
  list(simulation = kinship.data), "95%", "Synthetic kinship"
  )
assert.discrete.dot.style(
  kinship.plot, kinship$KINSHIP.CONTRAST.COLORS[kinship.contrasts],
  "Kinship contrasts"
  )

sfs <- helpers$sfs.R
sfs.contrasts <- c("AFR-ADX", "AFR-EUR", "ADX-EUR")
sfs.data <- make.discrete.contrast.data(sfs.contrasts, "minor.allele.count") %>%
  mutate(
    minor.allele.count = 1L,
    measure = factor("count", levels = c("count", "proportion"))
    )
sfs.plot <- sfs$make.sfs.contrast.plot(
  list(simulation = sfs.data), "95%", "Synthetic SFS"
  )
assert.discrete.dot.style(
  sfs.plot, sfs$SFS.CONTRAST.COLORS[sfs.contrasts], "SFS contrasts"
  )

# verify LD preserves ribbons and significant red-outlined points without legend.
ld <- helpers$ld_decay.R
ld.contrasts <- c("AFR-ADX", "AFR-EUR", "ADX-EUR")
ld.data <- tibble(
  distance_bin_bp = c(5000, 10000),
  contrast = factor("AFR-ADX", levels = ld.contrasts),
  difference = c(0.20, 0.15),
  ci.lower = c(0.05, 0.02),
  ci.upper = c(0.30, 0.25),
  point.color = "red"
  )
ld.plot <- ld$make.ld.simulation.contrast.plot(
  list(simulation = ld.data), "95", "population"
  )
assert.style(
  any(vapply(ld.plot$layers, function(layer) {
    inherits(layer$geom, "GeomRibbon")
    }, logical(1))),
  "LD contrasts do not retain confidence ribbons"
  )
assert.style(
  any(tolower(built.layer(ld.plot, "GeomPoint")$colour) %in%
    c("#ff0000", "red")),
  paste(
    "LD points do not retain significant red outlines:",
    paste(unique(built.layer(ld.plot, "GeomPoint")$colour), collapse = ", ")
    )
  )
assert.style(
  identical(ld.plot$theme$legend.position, "none"),
  "LD contrast facets retain a redundant legend"
  )


# render representative synthetic plots for visual inspection in a temporary path.
render.directory <- tempfile("plot-style-contract-")
dir.create(render.directory)
ggsave(file.path(render.directory, "ancestry.png"), ancestry.plot,
  width = 8, height = 4, dpi = 100)
ggsave(file.path(render.directory, "diversity.png"), diversity.plot,
  width = 8, height = 4, dpi = 100)
ggsave(file.path(render.directory, "kinship.png"), kinship.plot,
  width = 8, height = 4, dpi = 100)
ggsave(file.path(render.directory, "sfs.png"), sfs.plot,
  width = 8, height = 4, dpi = 100)
ggsave(file.path(render.directory, "ld.png"), ld.plot,
  width = 8, height = 4, dpi = 100)
assert.style(
  all(file.exists(file.path(render.directory,
    c("ancestry.png", "diversity.png", "kinship.png", "sfs.png", "ld.png")))),
  "Synthetic plot rendering did not create all representative images"
  )


cat("All synthetic plot-style contract tests passed\n")

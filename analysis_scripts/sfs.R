# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology
# Mooney Lab
# ___
# sfs.R
# ______________________________________________________________________________

# pattern: Mixed (unavoidable)
# Reason: This analysis script combines Parquet I/O, bootstrap calculations,
# and plot persistence.


# set up ----
library(tidyverse)
library(nanoparquet)
library(scales)


SIM.TC.DATA.DIR <- "~/scratch/OOA_NAAdmixture_2T12Consistent/stats"
SIMDOWN.TC.DATA.DIR <- "~/scratch/OOA_NAAdmixture_2T12ConsistentOnekgDownsample/stats"
SIM.LG.DATA.DIR <- "~/scratch/OOA_NAAdmixture_largeGrowth/stats"
SIMDOWN.LG.DATA.DIR <- "~/scratch/OOA_NAAdmixture_largeGrowthOnekgDownsample/stats"
EMPIRICAL.DATA.DIR <- "~/scratch/OOA_NAAdmixture_1kG/stats"
OUTPUT.DIR <- "/home1/karatas/proj/OOA_NAAdmixture_data"
CHROMOSOMES <- as.character(1:22)
SELECTED.CHROMOSOMES <- c("1")
DISPLAY.BIN.MAX <- 15
SFS.PROJECTION.ALLELE.COUNT <- 100
SOURCE.LEVELS <- c(
  "Simulation_2T12Consistent", "Simulation_2T12Consistent_simDown",
  "Simulation_largeGrowth", "Simulation_largeGrowth_simDown",
  "Empirical"
  )
SOURCE.LABELS <- setNames(
  c("T.C.", "T.C.D.", "L.G.", "L.G.D.", "Emp."), SOURCE.LEVELS
  )
POPULATION.LEVELS <- c("AFR", "ADX", "EUR", "YRI", "ASW", "CEU")
PLOT.CONFIGS <- list(
  tcd.1kg = SOURCE.LEVELS[c(2, 5)],
  tc.tcd.1kg = SOURCE.LEVELS[c(1, 2, 5)],
  all.datatypes.adx.asw = SOURCE.LEVELS
  )
SFS.FACET.LEVELS <- SELECTED.CHROMOSOMES
SFS.SERIES.PREFIXES <- c(
  Simulation_2T12Consistent = "TC",
  Simulation_2T12Consistent_simDown = "TC D.",
  Simulation_largeGrowth = "LG",
  Simulation_largeGrowth_simDown = "LG D.",
  Empirical = "empirical"
  )
SFS.SOURCE.POPULATIONS <- list(
  Simulation_2T12Consistent = POPULATION.LEVELS[1:3],
  Simulation_2T12Consistent_simDown = POPULATION.LEVELS[1:3],
  Simulation_largeGrowth = POPULATION.LEVELS[2],
  Simulation_largeGrowth_simDown = POPULATION.LEVELS[2],
  Empirical = POPULATION.LEVELS[4:6]
  )
SFS.SERIES.LEVELS <- unname(unlist(lapply(SOURCE.LEVELS, function(source) {
  paste(SFS.SERIES.PREFIXES[[source]], SFS.SOURCE.POPULATIONS[[source]])
  })))
SFS.COLORS <- c(
  "TC AFR" = "#9BD5F2",
  "TC ADX" = "#9A83CE",
  "TC EUR" = "#FBB4AE",
  "TC D. AFR" = "#56B4E9",
  "TC D. ADX" = "#6F55B5",
  "TC D. EUR" = "#FB8072",
  "LG ADX" = "#32146F",
  "LG D. ADX" = "#4B1FA8",
  "empirical YRI" = "#EEC4DC",
  "empirical ASW" = "#E44B8D",
  "empirical CEU" = "#BB437E"
  )
SFS.SERIES.LABELS <- c(
  "TC AFR" = "T.C. AFR",
  "TC ADX" = "T.C. ADX",
  "TC EUR" = "T.C. EUR",
  "TC D. AFR" = "T.C.D. AFR",
  "TC D. ADX" = "T.C.D. ADX",
  "TC D. EUR" = "T.C.D. EUR",
  "LG ADX" = "L.G. ADX",
  "LG D. ADX" = "L.G.D. ADX",
  "empirical YRI" = "YRI",
  "empirical ASW" = "ASW",
  "empirical CEU" = "CEU"
  )
SFS.BIN.WIDTH <- 1
SFS.DODGE <- position_dodge(width = SFS.BIN.WIDTH)
PLOT.BASE.SIZE <- 24
CATEGORICAL.BAR.LINEWIDTH <- 1
DENSE.BAR.WIDTH.MULTIPLIER <- 0.8
DENSE.BAR.LINEWIDTH <- 0.75
SFS.CONTRAST.BOOTSTRAP.REPLICATES <- 50L
SFS.CONTRAST.FAMILY.SIZE <- 45L
SFS.POPULATION.CONTRAST.SOURCES <- SOURCE.LEVELS[1:2]
SFS.POPULATION.CONTRASTS <- tribble(
  ~contrast, ~left.pop, ~right.pop,
  "AFR-ADX", "AFR", "ADX",
  "AFR-EUR", "AFR", "EUR",
  "ADX-EUR", "ADX", "EUR"
  )
SFS.SIMULATION.CONTRAST.SOURCES <- SOURCE.LEVELS[1:4]
SFS.SIMULATION.CONTRASTS <- tribble(
  ~contrast, ~left.source, ~right.source,
  "T.C. - T.C.D.", SOURCE.LEVELS[[1L]], SOURCE.LEVELS[[2L]],
  "L.G. - L.G.D.", SOURCE.LEVELS[[3L]], SOURCE.LEVELS[[4L]],
  "T.C. - L.G.", SOURCE.LEVELS[[1L]], SOURCE.LEVELS[[3L]]
  )
SFS.EMPIRICAL.CONTRASTS <- tribble(
  ~contrast, ~left.pop, ~right.pop,
  "YRI-ASW", "YRI", "ASW",
  "YRI-CEU", "YRI", "CEU",
  "ASW-CEU", "ASW", "CEU"
  )
SFS.CONTRAST.COLORS <- c(
  "AFR-ADX" = "#BFFBFF", "AFR-EUR" = "#16ACBD", "ADX-EUR" = "#00606F",
  "T.C. - T.C.D." = "#BFFBFF", "L.G. - L.G.D." = "#16ACBD",
  "T.C. - L.G." = "#00606F",
  "YRI-ASW" = "#BFFBFF", "YRI-CEU" = "#16ACBD", "ASW-CEU" = "#00606F"
  )
SFS.CONTRAST.LABELS <- c(
  "AFR-ADX" = "AFR - ADX", "AFR-EUR" = "AFR - EUR",
  "ADX-EUR" = "ADX - EUR",
  "T.C. - T.C.D." = "T.C. - T.C.D.",
  "L.G. - L.G.D." = "L.G. - L.G.D.", "T.C. - L.G." = "T.C. - L.G.",
  "YRI-ASW" = "YRI - ASW", "YRI-CEU" = "YRI - CEU",
  "ASW-CEU" = "ASW - CEU"
  )


# internal functions ----


# apply the common visual treatment to a completed plot.
apply.standard.plot.theme <- function(plot, legend.position = "top") {
  return(plot + theme_bw(base_size = PLOT.BASE.SIZE) + theme(
    legend.position = legend.position,
    legend.direction = "horizontal",
    legend.box = "horizontal",
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "plain")
    ))
  }


# return the canonical levels represented by a filtered plot view
order.active.levels <- function(values, canonical.levels) {
  active.levels <- canonical.levels[
    canonical.levels %in% as.character(values)
    ]
  return(active.levels)
  }


# summarize complete folded spectrum vectors with percentile intervals
summarize.bootstrap.interval <- function(
    data, grouping.columns, value.column
  ) {
  if (any(!is.finite(data[[value.column]]))) {
    stop("Simulation values must be finite")
    }
  summary <- data %>%
    group_by(across(all_of(setdiff(grouping.columns, "rep")))) %>%
    summarise(
      replicate.count = n_distinct(rep),
      mean = mean(.data[[value.column]]),
      lower = quantile(.data[[value.column]], 0.025, names = FALSE),
      upper = quantile(.data[[value.column]], 0.975, names = FALSE),
      .groups = "drop"
      )
  if (any(summary$replicate.count != 50L)) {
    stop("Simulation summaries require exactly 50 complete replicates")
    }
  return(summary)
  }


# retain source-specific populations before SFS normalization
apply.sfs.source.contract <- function(data) {
  retained <- data %>%
    filter(
      (data.type %in% c(
        "Simulation_2T12Consistent", "Simulation_2T12Consistent_simDown"
        ) & pop %in% c("AFR", "ADX", "EUR")) |
        (data.type %in% c(
          "Simulation_largeGrowth", "Simulation_largeGrowth_simDown"
          ) & pop == "ADX") |
        (data.type == "Empirical" & pop %in% c("YRI", "ASW", "CEU"))
      ) %>%
    mutate(data.type = factor(data.type, levels = SOURCE.LEVELS))
  return(retained)
  }


# require a table to contain its producer-owned SFS schema
check.sfs.columns <- function(data, required.columns, table.name) {
  missing.columns <- setdiff(required.columns, names(data))
  if (length(missing.columns) > 0L) {
    stop(paste(
      table.name,
      "is missing required columns:",
      paste(missing.columns, collapse = ", ")
      ))
    }
  return(invisible(NULL))
  }


# validate complete bins for every spectrum group
validate.complete.sfs <- function(data, allele.column, expected.bins) {
  invalid.count <- !is.finite(data$count) | data$count < 0
  if (any(invalid.count)) {
    stop("SFS counts must be finite and nonnegative")
    }
  bins <- data %>%
    group_by(data.set, rep, chrom, pop) %>%
    summarize(
      valid = identical(
        as.integer(sort(unique(.data[[allele.column]]))),
        as.integer(expected.bins)
        ),
      .groups = "drop"
      )
  if (!all(bins$valid)) {
    stop("SFS inputs must contain every expected allele-count bin")
    }
  return(data)
  }


# fold projected unfolded simulation spectra exactly once
fold.simulation.sfs <- function(data) {
  folded <- data %>%
    mutate(
      minor.allele.count = pmin(
        derived_allele_count,
        SFS.PROJECTION.ALLELE.COUNT - derived_allele_count
        )
      ) %>%
    group_by(data.set, data.type, rep, chrom, pop, minor.allele.count) %>%
    summarize(count = sum(count), .groups = "drop")
  return(folded)
  }


# add a genome spectrum by summing projected chromosome spectra
add.simulation.genome <- function(data) {
  genome <- data %>%
    group_by(data.set, data.type, rep, pop, derived_allele_count) %>%
    summarize(count = sum(count), .groups = "drop") %>%
    mutate(chrom = "all", .before = pop)
  return(bind_rows(data, genome))
  }


# add a genome spectrum by summing projected empirical chromosomes
add.empirical.genome <- function(data) {
  genome <- data %>%
    group_by(data.set, data.type, rep, pop, minor.allele.count) %>%
    summarize(count = sum(count), .groups = "drop") %>%
    mutate(chrom = "all", .before = pop)
  return(bind_rows(data, genome))
  }


# standardize producers, fold simulation once, and normalize full spectra
prepare.sfs.analysis <- function(simulation, simDown = NULL, empirical = NULL) {
  if (is.null(empirical)) {
    empirical <- simDown
    simDown <- NULL
    }
  if (!"data.set" %in% names(simulation)) {
    simulation$data.set <- ifelse(
      simulation$data.type %in% c(
        "Simulation_largeGrowth", "Simulation_largeGrowth_simDown"
        ),
      "largeGrowth", "2T12Consistent"
      )
    }
  check.sfs.columns(
    simulation,
    c("data.set", "rep", "chrom", "pop", "derived_allele_count", "count"),
    "Simulation SFS"
    )
  check.sfs.columns(
    empirical,
    c(
      "rep", "chrom", "pop", "minor_allele_count", "count",
      "projection_allele_count"
      ),
    "Empirical SFS"
    )
  if (!is.null(simDown)) {
    check.sfs.columns(
      simDown,
      c(
        "data.type", "rep", "chrom", "pop", "minor_allele_count",
        "count", "projection_allele_count"
        ),
      "simDown SFS"
      )
    invalid.projection <- is.na(simDown$projection_allele_count) |
      simDown$projection_allele_count != SFS.PROJECTION.ALLELE.COUNT
    if (any(invalid.projection)) {
      stop("simDown SFS projection metadata do not match the analysis")
      }
    }
  invalid.projection <- is.na(empirical$projection_allele_count) |
    empirical$projection_allele_count != SFS.PROJECTION.ALLELE.COUNT
  if (any(invalid.projection)) {
    stop("Empirical SFS projection metadata do not match the analysis")
    }

  simulation <- simulation %>%
    mutate(data.type = if_else(
      data.set == "2T12Consistent",
      "Simulation_2T12Consistent", "Simulation_largeGrowth"
      )) %>%
    mutate(chrom = as.character(chrom)) %>%
    validate.complete.sfs(
      "derived_allele_count",
      0:SFS.PROJECTION.ALLELE.COUNT
      ) %>%
    fold.simulation.sfs()
  if (!is.null(simDown)) {
    simDown <- simDown %>%
      mutate(
        data.type = as.character(data.type),
        data.set = if_else(
          grepl("largeGrowth", data.type),
          "largeGrowth", "2T12Consistent"
          ),
        chrom = as.character(chrom)
        ) %>%
      validate.complete.sfs(
        "minor_allele_count",
        0:(SFS.PROJECTION.ALLELE.COUNT / 2)
        ) %>%
      rename(minor.allele.count = minor_allele_count) %>%
      select(
        data.set, data.type, rep, chrom, pop, minor.allele.count, count
        )
    }
  empirical <- empirical %>%
    mutate(data.set = "empirical", data.type = "Empirical",
           chrom = as.character(chrom)) %>%
    validate.complete.sfs(
      "minor_allele_count",
      0:(SFS.PROJECTION.ALLELE.COUNT / 2)
      ) %>%
    rename(minor.allele.count = minor_allele_count) %>%
    add.empirical.genome() %>%
    select(
      data.set, data.type, rep, chrom, pop, minor.allele.count, count
      )

  prepared <- bind_rows(simulation, simDown, empirical) %>%
    apply.sfs.source.contract() %>%
    filter(chrom %in% c(SFS.FACET.LEVELS, "all")) %>%
    mutate(
      series = case_when(
        data.type == "Simulation_2T12Consistent" ~ paste("TC", pop),
        data.type == "Simulation_2T12Consistent_simDown" ~
          paste("TC D.", pop),
        data.type == "Simulation_largeGrowth" ~ paste("LG", pop),
        data.type == "Simulation_largeGrowth_simDown" ~
          paste("LG D.", pop),
        data.type == "Empirical" ~ paste("empirical", pop)
        ),
      pop = factor(pop, levels = POPULATION.LEVELS),
      chrom = factor(chrom, levels = c(SFS.FACET.LEVELS, "all")),
      series = factor(series, levels = SFS.SERIES.LEVELS)
      )
  if (any(is.na(prepared$series))) {
    stop("SFS series do not match the configured color keys")
    }

  prepared <- prepared %>%
    group_by(data.set, rep, chrom, pop, series) %>%
    mutate(
      segregating.total = sum(count[minor.allele.count %in% 1:50]),
      proportion = count / segregating.total
      ) %>%
    ungroup() %>%
    filter(minor.allele.count %in% 1:50) %>%
    select(-segregating.total)
  if (any(!is.finite(prepared$proportion))) {
    stop("Every SFS group must contain positive segregating-site mass")
    }
  prepared$data.type <- factor(prepared$data.type, levels = SOURCE.LEVELS)
  return(prepared)
  }


# summarize simulation replicates and retain empirical NA intervals
summarize.one.sfs.value <- function(data, value.column) {
  simulation <- data %>%
    filter(data.set != "empirical") %>%
    summarize.bootstrap.interval(
      c("data.set", "data.type", "rep", "chrom", "pop", "series",
        "minor.allele.count"),
      value.column
      )
  empirical <- data %>%
    filter(data.set == "empirical") %>%
    group_by(
      data.set, data.type, chrom, pop, series, minor.allele.count
      ) %>%
    summarize(
      mean = mean(.data[[value.column]]),
      .groups = "drop"
    ) %>%
    mutate(
      lower = NA_real_, upper = NA_real_, replicate.count = 1L
      )
  summary <- bind_rows(simulation, empirical) %>%
    mutate(
      chrom = factor(chrom, levels = SFS.FACET.LEVELS),
      series = factor(series, levels = SFS.SERIES.LEVELS)
      )
  return(summary)
  }


# return count and full-spectrum proportion summaries
summarize.sfs.analysis <- function(data) {
  summaries <- list(
    count = summarize.one.sfs.value(data, "count"),
    proportion = summarize.one.sfs.value(data, "proportion")
    )
  return(summaries)
  }


# require complete, unique, finite selected-chromosome contrast inputs
validate.sfs.contrast.input <- function(data, sources, populations) {
  required.columns <- c(
    "data.type", "rep", "chrom", "pop", "minor.allele.count",
    "count", "proportion"
    )
  check.sfs.columns(data, required.columns, "SFS contrast inputs")
  data <- data %>%
    filter(
      data.type %in% sources,
      pop %in% populations,
      as.character(chrom) %in% SELECTED.CHROMOSOMES,
      minor.allele.count %in% seq_len(DISPLAY.BIN.MAX)
      ) %>%
    mutate(
      data.type = as.character(data.type),
      chrom = as.character(chrom), pop = as.character(pop)
      )
  expected <- crossing(
    data.type = sources, rep = seq_len(50L),
    chrom = SELECTED.CHROMOSOMES, pop = populations,
    minor.allele.count = seq_len(DISPLAY.BIN.MAX)
    )
  counts <- data %>%
    count(data.type, rep, chrom, pop, minor.allele.count,
      name = "value.count")
  invalid <- expected %>%
    left_join(
      counts,
      by = c("data.type", "rep", "chrom", "pop", "minor.allele.count")
      ) %>%
    mutate(value.count = replace_na(value.count, 0L)) %>%
    filter(value.count != 1L)
  if (nrow(invalid) || any(!is.finite(data$count)) ||
      any(!is.finite(data$proportion))) {
    stop("SFS contrast inputs require exactly one finite value for every ",
      "source, replicate, population, and bin")
    }
  return(data)
  }


# bootstrap a paired left-minus-right selected-chromosome SFS difference
summarize.sfs.contrast.bootstrap <- function(
    left.values, right.values, bootstrap.replicates,
    family.size = SFS.CONTRAST.FAMILY.SIZE
  ) {
  if (length(left.values) != 50L || length(right.values) != 50L ||
      any(!is.finite(left.values)) || any(!is.finite(right.values))) {
    stop("SFS contrast bootstraps require 50 finite paired values")
    }
  draws <- replicate(bootstrap.replicates, {
    indices <- sample(seq_along(left.values), length(left.values),
      replace = TRUE)
    mean(left.values[indices] - right.values[indices])
    })
  nominal <- quantile(draws, c(0.025, 0.975), names = FALSE)
  bonferroni.quantile <- 0.05 / (2 * family.size)
  bonferroni <- quantile(
    draws, c(bonferroni.quantile, 1 - bonferroni.quantile), names = FALSE
    )
  return(tibble(
    difference = mean(left.values) - mean(right.values),
    ci.95.lower = nominal[[1L]], ci.95.upper = nominal[[2L]],
    bonferroni.ci.lower = bonferroni[[1L]],
    bonferroni.ci.upper = bonferroni[[2L]],
    bonferroni.quantile = bonferroni.quantile
    ))
  }


# construct simulation population contrast intervals for count and proportion
make.sfs.population.contrast.tables <- function(
    data, bootstrap.replicates = SFS.CONTRAST.BOOTSTRAP.REPLICATES,
    seed = 1L
  ) {
  data <- validate.sfs.contrast.input(
    data, SFS.POPULATION.CONTRAST.SOURCES, c("AFR", "ADX", "EUR")
    )
  set.seed(seed)
  tables <- map_dfr(SFS.POPULATION.CONTRAST.SOURCES, function(source) {
    map_dfr(seq_len(nrow(SFS.POPULATION.CONTRASTS)), function(index) {
      contrast <- SFS.POPULATION.CONTRASTS[index, ]
      map_dfr(c("count", "proportion"), function(measure) {
        map_dfr(seq_len(DISPLAY.BIN.MAX), function(bin) {
          left <- data %>% filter(
            data.type == source, pop == contrast$left.pop,
            minor.allele.count == bin
            ) %>% arrange(rep)
          right <- data %>% filter(
            data.type == source, pop == contrast$right.pop,
            minor.allele.count == bin
            ) %>% arrange(rep)
          if (!identical(left$rep, right$rep)) {
            stop("Paired SFS contrast replicate IDs must match")
            }
          bind_cols(
            tibble(
              data.type = source, minor.allele.count = bin,
              contrast = contrast$contrast, measure = measure
              ),
            summarize.sfs.contrast.bootstrap(
              left[[measure]], right[[measure]], bootstrap.replicates
              )
            )
          })
        })
      })
    })
  return(tables)
  }


# construct paired ADX source contrast intervals for count and proportion
make.sfs.simulation.contrast.tables <- function(
    data, bootstrap.replicates = SFS.CONTRAST.BOOTSTRAP.REPLICATES,
    seed = 1L
  ) {
  data <- validate.sfs.contrast.input(
    data, SFS.SIMULATION.CONTRAST.SOURCES, "ADX"
    )
  set.seed(seed)
  tables <- map_dfr(seq_len(nrow(SFS.SIMULATION.CONTRASTS)), function(index) {
    contrast <- SFS.SIMULATION.CONTRASTS[index, ]
    map_dfr(c("count", "proportion"), function(measure) {
      map_dfr(seq_len(DISPLAY.BIN.MAX), function(bin) {
        left <- data %>% filter(
          data.type == contrast$left.source, minor.allele.count == bin
          ) %>% arrange(rep)
        right <- data %>% filter(
          data.type == contrast$right.source, minor.allele.count == bin
          ) %>% arrange(rep)
        if (!identical(left$rep, right$rep)) {
          stop("Paired SFS contrast replicate IDs must match")
          }
        bind_cols(
          tibble(
            minor.allele.count = bin, contrast = contrast$contrast,
            measure = measure
            ),
          summarize.sfs.contrast.bootstrap(
            left[[measure]], right[[measure]], bootstrap.replicates
            )
          )
        })
      })
    })
  return(tables)
  }


# calculate direct empirical selected-chromosome population differences
make.sfs.empirical.contrast.table <- function(data) {
  data <- data %>%
    filter(
      data.type == "Empirical", as.character(chrom) %in%
        SELECTED.CHROMOSOMES,
      pop %in% c("YRI", "ASW", "CEU"),
      minor.allele.count %in% seq_len(DISPLAY.BIN.MAX)
      )
  expected <- crossing(
    pop = c("YRI", "ASW", "CEU"),
    minor.allele.count = seq_len(DISPLAY.BIN.MAX)
    )
  counts <- data %>% count(pop, minor.allele.count, name = "value.count")
  invalid <- expected %>% left_join(
    counts, by = c("pop", "minor.allele.count")
    ) %>% mutate(value.count = replace_na(value.count, 0L)) %>%
    filter(value.count != 1L)
  if (nrow(invalid) || any(!is.finite(data$count)) ||
      any(!is.finite(data$proportion))) {
    stop("Empirical SFS contrast inputs require exactly one finite value")
    }
  tables <- map_dfr(seq_len(nrow(SFS.EMPIRICAL.CONTRASTS)), function(index) {
    contrast <- SFS.EMPIRICAL.CONTRASTS[index, ]
    contrast.label <- contrast$contrast
    map_dfr(c("count", "proportion"), function(measure) {
      left <- data %>% filter(pop == contrast$left.pop) %>%
        select(minor.allele.count, left.value = all_of(measure))
      right <- data %>% filter(pop == contrast$right.pop) %>%
        select(minor.allele.count, right.value = all_of(measure))
      inner_join(left, right, by = "minor.allele.count") %>%
        transmute(
          minor.allele.count, contrast = contrast.label,
          populations = contrast.label, measure = measure,
          difference = left.value - right.value
          )
      })
    })
  return(tables)
  }


# select one interval family and establish contrast plotting order
prepare.sfs.contrast.plot.data <- function(
    data, interval.type = c("95", "bonferroni"), source.group = NULL
  ) {
  interval.type <- match.arg(interval.type)
  lower.column <- if (interval.type == "95") {
    "ci.95.lower"
    } else {
    "bonferroni.ci.lower"
    }
  upper.column <- if (interval.type == "95") {
    "ci.95.upper"
    } else {
    "bonferroni.ci.upper"
    }
  if (identical(source.group, "tc.tcd")) {
    source.group <- SFS.POPULATION.CONTRAST.SOURCES
    }
  contrasts <- if ("data.type" %in% names(data)) {
    SFS.POPULATION.CONTRASTS$contrast
    } else {
    SFS.SIMULATION.CONTRASTS$contrast
    }
  displayed <- if (is.null(source.group)) {
    data
    } else {
    data %>% filter(data.type %in% source.group)
    }
  displayed <- displayed %>%
    transmute(
      across(any_of("data.type")), minor.allele.count, contrast, measure,
      difference, ci.lower = .data[[lower.column]],
      ci.upper = .data[[upper.column]]
      ) %>%
    mutate(
      contrast = factor(contrast, levels = contrasts),
      measure = factor(measure, levels = c("count", "proportion")),
      outside.zero = ci.lower > 0 | ci.upper < 0,
      point.color = if_else(
        outside.zero, "red", "black"
        )
      )
  if ("data.type" %in% names(displayed)) {
    displayed <- displayed %>% mutate(data.type = factor(
      data.type, levels = SFS.POPULATION.CONTRAST.SOURCES
      ))
    }
  return(list(
    simulation = displayed,
    red.markers = filter(displayed, outside.zero)
    ))
  }


# build selected-chromosome contrast plot with a zero reference and dodging
make.sfs.contrast.plot <- function(data, interval.label, title) {
  simulation <- data$simulation
  contrasts <- levels(simulation$contrast)
  active.contrasts <- contrasts[
    contrasts %in% unique(as.character(simulation$contrast))
    ]
  plot <- ggplot(
    simulation,
    aes(minor.allele.count, difference, group = contrast)
    ) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1) +
    geom_errorbar(
      aes(ymin = ci.lower, ymax = ci.upper), position = SFS.DODGE, width = 0,
      color = "black", linewidth = CATEGORICAL.BAR.LINEWIDTH,
      show.legend = FALSE
      ) +
    geom_col(
      aes(fill = contrast, color = point.color), position = SFS.DODGE,
      width = SFS.BIN.WIDTH * DENSE.BAR.WIDTH.MULTIPLIER,
      linewidth = CATEGORICAL.BAR.LINEWIDTH
      ) +
    scale_x_continuous(breaks = seq_len(DISPLAY.BIN.MAX)) +
    scale_color_manual(values = c(black = "black", red = "red"),
      guide = "none") +
    scale_fill_manual(
      values = SFS.CONTRAST.COLORS[active.contrasts],
      breaks = active.contrasts,
      labels = SFS.CONTRAST.LABELS[active.contrasts]
      ) +
    labs(
      x = "Minor allele count bin", y = "Left population/source minus right",
      color = NULL, title = title, subtitle = paste(interval.label, "interval")
      ) +
    guides(
      fill = guide_legend(
        nrow = 1, byrow = TRUE,
        override.aes = list(
          shape = 22,
          fill = unname(SFS.CONTRAST.COLORS[active.contrasts]),
          color = "black"
          )
        )
      ) +
    theme()
  if ("data.type" %in% names(simulation)) {
    plot <- plot + facet_grid(
      rows = vars(data.type, measure), scales = "free_y",
      labeller = labeller(data.type = SOURCE.LABELS)
      )
    } else {
    plot <- plot + facet_grid(measure ~ ., scales = "free_y")
    }
  return(apply.standard.plot.theme(plot))
  }


# build direct empirical selected-chromosome contrast reference plot
make.sfs.empirical.contrast.plot <- function(data) {
  contrasts <- SFS.EMPIRICAL.CONTRASTS$contrast
  active.contrasts <- contrasts[
    contrasts %in% unique(as.character(data$contrast))
    ]
  plot <- ggplot(
    data %>% mutate(
      contrast = factor(contrast, levels = contrasts),
      measure = factor(measure, levels = c("count", "proportion"))
      ),
    aes(minor.allele.count, difference, group = contrast)
    ) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1) +
    geom_col(
      aes(minor.allele.count, difference, fill = contrast, group = contrast),
      position = SFS.DODGE,
      width = SFS.BIN.WIDTH * DENSE.BAR.WIDTH.MULTIPLIER,
      color = "black", linewidth = CATEGORICAL.BAR.LINEWIDTH,
      inherit.aes = FALSE
      ) +
    scale_x_continuous(breaks = seq_len(DISPLAY.BIN.MAX)) +
    scale_fill_manual(
      values = SFS.CONTRAST.COLORS[active.contrasts],
      breaks = active.contrasts,
      labels = SFS.CONTRAST.LABELS[active.contrasts]
      ) +
    labs(
      x = "Minor allele count bin", y = "Left population minus right",
      fill = NULL, title = "Empirical population contrasts"
      ) +
    guides(
      fill = guide_legend(
        nrow = 1, byrow = TRUE,
        override.aes = list(
          shape = 22,
          fill = unname(SFS.CONTRAST.COLORS[active.contrasts]),
          color = "black"
          )
        )
      ) +
    facet_grid(measure ~ ., scales = "free_y") +
    theme()
  return(apply.standard.plot.theme(plot))
  }


# retain one configured SFS view
filter.plot.view <- function(data, data.types, tag) {
  if (!tag %in% names(PLOT.CONFIGS)) {
    stop("Unsupported SFS plot tag: ", tag)
    }
  if (!identical(data.types, PLOT.CONFIGS[[tag]])) {
    stop("SFS data types do not match the configured tag")
    }
  filtered <- data %>%
    filter(as.character(data.type) %in% data.types) %>%
    filter(
      tag != "all.datatypes.adx.asw" |
        (data.type != "Empirical" & pop == "ADX") |
        (data.type == "Empirical" & pop == "ASW")
      )
  active.sources <- order.active.levels(filtered$data.type, SOURCE.LEVELS)
  active.populations <- order.active.levels(
    filtered$pop, POPULATION.LEVELS
    )
  filtered <- filtered %>%
    mutate(
      data.type = factor(as.character(data.type), levels = active.sources),
      pop = factor(as.character(pop), levels = active.populations),
      series = factor(
        as.character(series),
        levels = order.active.levels(series, SFS.SERIES.LEVELS)
        ),
      chrom = factor(
        as.character(chrom),
        levels = c(
          SFS.FACET.LEVELS[SFS.FACET.LEVELS %in% chrom],
          if ("all" %in% chrom) "all"
          )
        )
      ) %>%
    arrange(data.type, series) %>%
    droplevels()
  return(filtered)
  }


# build one scoped shared dodged-bar SFS plot
make.sfs.plot <- function(
    data, value.column, y.label, pseudo.log, data.types, tag,
    show.all = FALSE
  ) {
  displayed <- filter.plot.view(data, data.types, tag) %>%
    filter(
      as.character(chrom) %in% SELECTED.CHROMOSOMES |
        (show.all & data.type == "Empirical" & chrom == "all")
      ) %>%
    filter(minor.allele.count <= DISPLAY.BIN.MAX)
  series.keys <- levels(displayed$series)
  plot <- ggplot(
    displayed,
    aes(
      x = minor.allele.count,
      y = .data[[value.column]],
      fill = series,
      group = series
      )
    ) +
    geom_col(
      position = SFS.DODGE,
      width = SFS.BIN.WIDTH * DENSE.BAR.WIDTH.MULTIPLIER,
      color = "black", linewidth = DENSE.BAR.LINEWIDTH
      ) +
    geom_errorbar(
      aes(ymin = lower, ymax = upper),
      position = SFS.DODGE,
      width = 0, linewidth = DENSE.BAR.LINEWIDTH,
      na.rm = TRUE
      ) +
    scale_x_continuous(
      breaks = seq_len(DISPLAY.BIN.MAX),
      limits = c(0.5, DISPLAY.BIN.MAX + 0.5)
      ) +
    scale_fill_manual(
      values = SFS.COLORS[series.keys],
      breaks = series.keys,
      labels = SFS.SERIES.LABELS[series.keys],
      limits = series.keys,
      drop = TRUE
      ) +
    labs(
      x = "Minor allele count",
      y = y.label,
      title = paste(
        if (y.label == "Projected site count") "SFS counts" else "SFS",
        tag, sep = ": "
        ),
      subtitle = paste0(
        "Projected to ", SFS.PROJECTION.ALLELE.COUNT, " alleles"
        ),
      fill = NULL
      ) +
    guides(fill = guide_legend(order = 1, nrow = 1, byrow = TRUE))
  if (pseudo.log) {
    plot <- plot +
      scale_y_log10(
        breaks = c(1e4, 2.5e4, 5e4, 1e5, 2.5e5, 5e5, 1e6),
        labels = label_number(
          scale_cut = cut_short_scale()
          )
        )
    } else {
    plot <- plot + scale_y_continuous()
    }
  return(apply.standard.plot.theme(plot))
  }


# read every available projected chromosome file from one source
read.sfs.chromosomes <- function(directory, chromosomes, file.family.label) {
  paths <- file.path(
    path.expand(directory),
    paste0("sfs.chr", chromosomes, ".parquet")
    )
  missing <- !file.exists(paths)
  if (any(missing)) {
    warning(
      paste0(
        file.family.label, " is unavailable for chromosomes: ",
        paste(chromosomes[missing], collapse = ", ")
        ),
      call. = FALSE
      )
    }
  paths <- paths[!missing]
  chromosomes <- chromosomes[!missing]
  read.one <- function(path, chrom) {
    table <- read_parquet(path)
    table$chrom <- chrom
    return(table)
    }
  if (requireNamespace("furrr", quietly = TRUE)) {
    previous.plan <- future::plan()
    on.exit(future::plan(previous.plan), add = TRUE)
    future::plan(future::multisession)
    tables <- furrr::future_map2_dfr(
      paths, chromosomes, read.one,
      .options = furrr::furrr_options(seed = TRUE)
      )
    } else {
    tables <- purrr::map2_dfr(paths, chromosomes, read.one)
    }
  return(tables)
  }


# read standardized projected chromosome producer outputs
read.sfs.inputs <- function(
    sim.tc.dir, simDown.tc.dir, sim.lg.dir,
    simDown.lg.dir, empirical.dir, chromosomes
  ) {
  simulation <- bind_rows(
    read.sfs.chromosomes(
      sim.tc.dir, chromosomes, "TC SFS files"
      ) %>%
      mutate(data.set = "2T12Consistent"),
    read.sfs.chromosomes(
      sim.lg.dir, chromosomes, "LG SFS files"
      ) %>%
      mutate(data.set = "largeGrowth")
    )
  simDown <- bind_rows(
    read.sfs.chromosomes(
      simDown.tc.dir, chromosomes, "TC D. SFS files"
      ) %>%
      mutate(data.type = "Simulation_2T12Consistent_simDown"),
    read.sfs.chromosomes(
      simDown.lg.dir, chromosomes, "LG D. SFS files"
      ) %>%
      mutate(data.type = "Simulation_largeGrowth_simDown")
    )
  empirical <- read.sfs.chromosomes(
    empirical.dir, chromosomes, "Empirical SFS files"
    )
  return(list(
    simulation = simulation,
    simDown = simDown,
    empirical = empirical
    ))
  }


# analysis ----


sfs.inputs <- read.sfs.inputs(
  SIM.TC.DATA.DIR,
  SIMDOWN.TC.DATA.DIR,
  SIM.LG.DATA.DIR,
  SIMDOWN.LG.DATA.DIR,
  EMPIRICAL.DATA.DIR,
  CHROMOSOMES
  )
sfs.data <- prepare.sfs.analysis(
  sfs.inputs$simulation,
  sfs.inputs$simDown,
  sfs.inputs$empirical
  )
sfs.summaries <- summarize.sfs.analysis(sfs.data)
sfs.population.contrast.tables <- make.sfs.population.contrast.tables(sfs.data)
sfs.simulation.contrast.tables <- make.sfs.simulation.contrast.tables(sfs.data)
sfs.empirical.contrast.table <- make.sfs.empirical.contrast.table(sfs.data)
sfs.population.contrast.plots <- list(
  tc.tcd = list(
    `95` = make.sfs.contrast.plot(
      prepare.sfs.contrast.plot.data(
        sfs.population.contrast.tables, "95",
        SFS.POPULATION.CONTRAST.SOURCES
        ),
      "95%", "T.C. and T.C.D. population contrasts"
      ),
    bonferroni = make.sfs.contrast.plot(
      prepare.sfs.contrast.plot.data(
        sfs.population.contrast.tables, "bonferroni",
        SFS.POPULATION.CONTRAST.SOURCES
        ),
      "Bonferroni", "T.C. and T.C.D. population contrasts"
      )
    ),
  empirical = make.sfs.empirical.contrast.plot(sfs.empirical.contrast.table)
  )
sfs.simulation.contrast.plots <- list(
  `95` = make.sfs.contrast.plot(
    prepare.sfs.contrast.plot.data(sfs.simulation.contrast.tables, "95"),
    "95%", "ADX simulation source contrasts"
    ),
  bonferroni = make.sfs.contrast.plot(
    prepare.sfs.contrast.plot.data(
      sfs.simulation.contrast.tables, "bonferroni"
      ),
    "Bonferroni", "ADX simulation source contrasts"
    )
  )

sfs.bootstrap.count.plots <- imap(list(
  tcd.1kg = PLOT.CONFIGS$tcd.1kg,
  all.datatypes.adx.asw = PLOT.CONFIGS$all.datatypes.adx.asw
  ), function(data.types, tag) {
  return(make.sfs.plot(
    sfs.summaries$count, "mean", "Projected site count", TRUE,
    data.types, tag,
    show.all = FALSE
    ))
  })
sfs.bootstrap.proportion.plots <- imap(list(
  tcd.1kg = PLOT.CONFIGS$tcd.1kg,
  all.datatypes.adx.asw = PLOT.CONFIGS$all.datatypes.adx.asw
  ), function(data.types, tag) {
  return(make.sfs.plot(
    sfs.summaries$proportion, "mean",
    "Proportion of segregating sites", FALSE, data.types,
    tag,
    show.all = FALSE
    ))
  })

# persist chromosome-1 bootstrap count and proportion plots before printing.
dir.create(OUTPUT.DIR, recursive = TRUE, showWarnings = FALSE)
saveRDS(sfs.bootstrap.count.plots$tcd.1kg, file.path(
  OUTPUT.DIR, "sfs.bootstrap.count.tcd.1kg.rds"
  ))
saveRDS(sfs.bootstrap.count.plots$all.datatypes.adx.asw, file.path(
  OUTPUT.DIR, "sfs.bootstrap.count.all.datatypes.adx.asw.rds"
  ))
saveRDS(sfs.bootstrap.proportion.plots$tcd.1kg, file.path(
  OUTPUT.DIR, "sfs.bootstrap.proportion.tcd.1kg.rds"
  ))
saveRDS(sfs.bootstrap.proportion.plots$all.datatypes.adx.asw, file.path(
  OUTPUT.DIR, "sfs.bootstrap.proportion.all.datatypes.adx.asw.rds"
  ))
saveRDS(sfs.population.contrast.plots$tc.tcd$`95`, file.path(
  OUTPUT.DIR, "sfs.population.contrasts.tc.tcd.95.rds"
  ))
saveRDS(sfs.population.contrast.plots$tc.tcd$bonferroni, file.path(
  OUTPUT.DIR, "sfs.population.contrasts.tc.tcd.bonferroni.rds"
  ))
saveRDS(sfs.population.contrast.plots$empirical, file.path(
  OUTPUT.DIR, "sfs.population.contrasts.empirical.rds"
  ))
saveRDS(sfs.simulation.contrast.plots$`95`, file.path(
  OUTPUT.DIR, "sfs.simulation.contrasts.95.rds"
  ))
saveRDS(sfs.simulation.contrast.plots$bonferroni, file.path(
  OUTPUT.DIR, "sfs.simulation.contrasts.bonferroni.rds"
  ))

print(sfs.bootstrap.count.plots$tcd.1kg)
print(sfs.bootstrap.count.plots$all.datatypes.adx.asw)
print(sfs.bootstrap.proportion.plots$tcd.1kg)
print(sfs.bootstrap.proportion.plots$all.datatypes.adx.asw)
# print(sfs.population.contrast.plots$tc.tcd$`95`)
print(sfs.population.contrast.plots$tc.tcd$bonferroni)
print(sfs.population.contrast.plots$empirical)
# print(sfs.simulation.contrast.plots$`95`)
print(sfs.simulation.contrast.plots$bonferroni)

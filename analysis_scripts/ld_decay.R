# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology
# Mooney Lab
# ___
# ld_decay.R
# ______________________________________________________________________________


# set up ----
options(scipen = 999)
library(tidyverse)
library(glue)
library(nanoparquet)


SIM.TC.DATA.DIR <- "~/scratch/OOA_NAAdmixture_2T12Consistent/stats"
SIMDOWN.TC.DATA.DIR <- "~/scratch/OOA_NAAdmixture_2T12ConsistentOnekgDownsample/stats"
SIM.LG.DATA.DIR <- "~/scratch/OOA_NAAdmixture_largeGrowth/stats"
SIMDOWN.LG.DATA.DIR <- "~/scratch/OOA_NAAdmixture_largeGrowthOnekgDownsample/stats"
EMPIRICAL.DATA.DIR <- "~/scratch/OOA_NAAdmixture_1kG/stats"
OUTPUT.DIR <- "/home1/karatas/proj/OOA_NAAdmixture_data"
CHROMOSOMES <- as.character(1:22)
SELECTED.CHROMOSOMES <- c("1")
SOURCE.LEVELS <- c(
  "Simulation_2T12Consistent", "Simulation_2T12Consistent_simDown",
  "Simulation_largeGrowth", "Simulation_largeGrowth_simDown",
  "Empirical"
  )
SOURCE.DISPLAY.LEVELS <- c("T.C.", "T.C.D.", "L.G.", "L.G.D.", "Emp.")
SOURCE.LABELS <- setNames(SOURCE.DISPLAY.LEVELS, SOURCE.LEVELS)
POPULATION.LEVELS <- c("AFR", "ADX", "EUR", "YRI", "ASW", "CEU")
PLOT.CONFIGS <- list(
  tcd.1kg = SOURCE.LEVELS[c(2, 5)],
  tc.tcd.1kg = SOURCE.LEVELS[c(1, 2, 5)],
  all.datatypes.adx.asw = SOURCE.LEVELS
  )
BOOTSTRAP.LEGEND.VIEWS <- c("all.lines", "role.interval")
PLOT.BASE.SIZE <- 24
LD.CONTRAST.LINEWIDTH <- 1
LD.PLOT.DISTANCE.BINS <- seq(5000, 200000, by = 5000)
LD.POPULATION.CONTRAST.BOOTSTRAP.REPLICATES <- 1000L
# per displayed comparison-set Bonferroni family: 3 contrasts x 40 bins.
LD.POPULATION.CONTRAST.FAMILY.SIZE <- 120L
LD.POPULATION.CONTRAST.SOURCES <- SOURCE.LEVELS[c(1L, 2L)]
LD.SIMULATION.CONTRAST.SOURCES <- SOURCE.LEVELS[seq_len(4L)]
LD.POPULATION.CONTRASTS <- tribble(
  ~contrast, ~simulation.left, ~simulation.right,
  ~empirical.left, ~empirical.right,
  "AFR-ADX", "AFR", "ADX", "YRI", "ASW",
  "AFR-EUR", "AFR", "EUR", "YRI", "CEU",
  "ADX-EUR", "ADX", "EUR", "ASW", "CEU"
  )
LD.SIMULATION.CONTRASTS <- tribble(
  ~contrast, ~left.source, ~right.source, ~paired,
  "T.C. - T.C.D.", SOURCE.LEVELS[[1L]], SOURCE.LEVELS[[2L]], TRUE,
  "L.G. - L.G.D.", SOURCE.LEVELS[[3L]], SOURCE.LEVELS[[4L]], TRUE,
  "T.C. - L.G.", SOURCE.LEVELS[[1L]], SOURCE.LEVELS[[3L]], FALSE
  )
PLOT.STYLES <- list(
  population.colors = c(
    AFR = "#56B4E9", ADX = "#4B1FA8", EUR = "#fb8072",
    YRI = "#eec4dc", ASW = "#e44b8d", CEU = "#bb437e"
    ),
  source.colors = c(
    Simulation_2T12Consistent = "#9A83CE",
    Simulation_2T12Consistent_simDown = "#6F55B5",
    Simulation_largeGrowth = "#32146F",
    Simulation_largeGrowth_simDown = "#4B1FA8"
    ),
  contrast.colors = c(
    `AFR-ADX` = "#BFFBFF", `AFR-EUR` = "#16ACBD",
    `ADX-EUR` = "#00606F"
    ),
  simulation.contrast.colors = c(
    "T.C. - T.C.D." = "#BFFBFF", "L.G. - L.G.D." = "#16ACBD",
    "T.C. - L.G." = "#00606F"
    ),
  series.labels = SOURCE.LABELS
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


# summarize complete pooled replicate curves with percentile intervals
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


# retain source-specific populations before pooling sufficient statistics
apply.ld.source.contract <- function(data) {
  retained <- data %>%
    filter(
      (data.type %in% c(
        "Simulation_2T12Consistent", "Simulation_2T12Consistent_simDown"
        ) & pop %in% c("AFR", "ADX", "EUR")) |
        (data.type %in% c(
          "Simulation_largeGrowth", "Simulation_largeGrowth_simDown"
          ) & pop == "ADX") |
        (data.type == "Empirical" & pop %in% c("AFR", "ADX", "EUR",
                                                "YRI", "ASW", "CEU"))
      ) %>%
    mutate(data.type = factor(data.type, levels = SOURCE.LEVELS))

  return(retained)
  }


# map source population labels to shared population roles
add.population.roles <- function(data) {
  data <- data %>%
    mutate(role = case_when(
      pop %in% c("AFR", "YRI") ~ "AFR",
      pop %in% c("ADX", "ASW") ~ "ADX",
      pop %in% c("EUR", "CEU") ~ "EUR",
      TRUE ~ NA_character_
      ))
  if (any(is.na(data$role))) {
    stop("LD data contain an unsupported population label")
    }

  return(data)
  }


# standardize one LD table to the shared analysis schema
normalize.ld.table <- function(data, data.type.input) {
  data <- data %>%
    mutate(
      rep = as.numeric(rep),
      chrom = as.character(chrom),
      pop = as.character(pop),
      distance_bin_bp = as.numeric(distance_bin_bp),
      mean_r2 = as.numeric(mean_r2),
      sum_r2 = as.numeric(sum_r2),
      n_pairs = as.numeric(n_pairs),
      data.type = data.type.input
      ) %>%
    apply.ld.source.contract() %>%
    add.population.roles()

  return(data)
  }


# read chromosome-labelled LD files for one source
read.ld.chromosomes <- function(
    data.directory, chromosomes, data.type.input
  ) {
  paths <- file.path(
    path.expand(data.directory),
    glue::glue("ld_decay.chr{chromosomes}.parquet")
    )
  missing <- !file.exists(paths)
  if (any(missing)) {
    warning(
      paste0(
        data.type.input,
        " LD files are unavailable for chromosomes: ",
        paste(chromosomes[missing], collapse = ", ")
        ),
      call. = FALSE
      )
    }
  paths <- paths[!missing]
  chromosomes <- chromosomes[!missing]
  data <- map2_dfr(paths, chromosomes, function(path, chrom) {
    table <- read_parquet(path)
    table$chrom <- chrom
    return(normalize.ld.table(table, data.type.input))
    })

  return(data)
  }


# read the empirical genome-wide LD table
read.empirical.ld.genome <- function(data.directory) {
  data <- read_parquet(file.path(
    path.expand(data.directory), "ld_decay.parquet"
    ))
  data$chrom <- "all"
  data <- normalize.ld.table(data, "Empirical")

  return(data)
  }


# reconstruct chromosome or genome curves from pooled sums and pair counts
pool.ld.curves <- function(data, include.chromosome) {
  if (!"role" %in% names(data)) data <- add.population.roles(data)
  grouping.columns <- c(
    "data.type", "rep", "pop", "role", "distance_bin_bp"
    )
  if (include.chromosome) {
    grouping.columns <- c(grouping.columns, "chrom")
    }
  pooled <- data %>%
    group_by(across(all_of(grouping.columns))) %>%
    summarise(
      sum.r2 = sum(sum_r2),
      n.pairs = sum(n_pairs),
      chromosome.count = if (include.chromosome) {
        NA_integer_
        } else {
        n_distinct(chrom)
        },
      .groups = "drop"
      ) %>%
    mutate(
      mean.r2 = if_else(n.pairs > 0, sum.r2 / n.pairs, NA_real_)
      )
  if (!include.chromosome) pooled$chrom <- "all"

  return(pooled)
  }


# summarize replicate simulation curves and fixed empirical curves
summarize.ld.curves <- function(data, chromosomes) {
  chromosomes <- as.character(chromosomes)
  if (!length(chromosomes) || any(!chromosomes %in% CHROMOSOMES)) {
    stop("LD summarization requires autosomal chromosomes")
    }
  scoped <- data %>% filter(as.character(chrom) %in% chromosomes)
  if (!nrow(scoped)) {
    stop("LD data do not contain any requested chromosomes")
    }
  if (any(!as.character(scoped$chrom) %in% chromosomes)) {
    stop("LD summary contains data outside the requested chromosomes")
    }
  simulation <- scoped %>%
    filter(data.type != "Empirical") %>%
    summarize.bootstrap.interval(
      c("data.type", "rep", "pop", "role", "chrom", "distance_bin_bp"),
      "mean.r2"
      )
  empirical <- scoped %>%
    filter(data.type == "Empirical") %>%
    transmute(
      data.type, pop, role, chrom, distance_bin_bp,
      mean = mean.r2, lower = NA_real_, upper = NA_real_, replicate.count = 1L,
      chromosome.count
      ) %>%
    distinct()
  empirical.genome <- data %>%
    filter(data.type == "Empirical", as.character(chrom) == "all") %>%
    transmute(
      data.type, pop, role, chrom, distance_bin_bp,
      mean = mean.r2, lower = NA_real_, upper = NA_real_, replicate.count = 1L,
      chromosome.count
      ) %>%
    distinct()
  summary <- bind_rows(simulation, empirical, empirical.genome) %>%
    mutate(
      role = factor(role, levels = POPULATION.LEVELS[1:3]),
      chrom = factor(as.character(chrom), levels = c(chromosomes, "all")),
      data.type = factor(
        data.type,
        levels = SOURCE.LEVELS
        )
      ) %>%
    filter(!is.na(mean))

  return(summary)
  }


# validate complete selected-chromosome simulation LD curves
validate.ld.contrast.simulation <- function(data, sources, populations) {
  required.columns <- c(
    "data.type", "rep", "chrom", "pop", "distance_bin_bp", "mean.r2"
    )
  missing.columns <- setdiff(required.columns, names(data))
  if (length(missing.columns)) {
    stop("LD contrast data are missing columns: ",
      paste(missing.columns, collapse = ", "))
    }
  data <- data %>%
    filter(
      data.type %in% sources,
      as.character(chrom) %in% SELECTED.CHROMOSOMES,
      pop %in% populations,
      distance_bin_bp %in% LD.PLOT.DISTANCE.BINS
      ) %>%
    mutate(
      data.type = as.character(data.type), chrom = as.character(chrom),
      pop = as.character(pop), distance_bin_bp = as.numeric(distance_bin_bp)
      )
  if (any(!is.finite(data$mean.r2))) {
    stop("LD contrast simulation values must be finite")
    }
  expected <- crossing(
    data.type = sources, rep = seq_len(50L), chrom = SELECTED.CHROMOSOMES,
    pop = populations, distance_bin_bp = LD.PLOT.DISTANCE.BINS
    )
  counts <- expected %>%
    left_join(
      data %>%
        count(data.type, rep, chrom, pop, distance_bin_bp,
          name = "value.count"),
      by = c("data.type", "rep", "chrom", "pop", "distance_bin_bp")
      ) %>%
    mutate(value.count = replace_na(value.count, 0L)) %>%
    filter(value.count != 1L)
  unexpected <- anti_join(
    data, expected,
    by = c("data.type", "rep", "chrom", "pop", "distance_bin_bp")
    )
  if (nrow(counts) || nrow(unexpected)) {
    stop("LD population contrast inputs require exactly one finite value for ",
      "every source, replicate, selected chromosome, population, and ",
      "displayed distance bin")
    }
  return(data)
  }


# validate paired replicate IDs for one set of populations or sources
validate.ld.paired.ids <- function(data, grouping.columns, value.column,
    error.message) {
  ids <- data %>%
    group_by(across(all_of(c(grouping.columns, value.column)))) %>%
    summarise(rep.ids = list(sort(rep)), .groups = "drop") %>%
    pivot_wider(names_from = all_of(value.column), values_from = rep.ids)
  value.names <- setdiff(names(ids), grouping.columns)
  for (row in seq_len(nrow(ids))) {
    values <- ids[row, value.names]
    if (any(vapply(values, is.null, logical(1))) ||
        !all(vapply(values, identical, logical(1), values[[1L]]))) {
      stop(error.message)
      }
    }
  return(invisible(NULL))
  }


# bootstrap one population difference from paired replicate curves
summarize.ld.population.contrast.bootstrap <- function(
    left.values, right.values, bootstrap.replicates,
    family.size = LD.POPULATION.CONTRAST.FAMILY.SIZE
  ) {
  if (length(left.values) != 50L || length(right.values) != 50L ||
      any(!is.finite(left.values)) || any(!is.finite(right.values))) {
    stop("Paired LD population contrast bootstraps require 50 finite values")
    }
  draws <- replicate(bootstrap.replicates, {
    indices <- sample(seq_along(left.values), length(left.values),
      replace = TRUE)
    mean(left.values[indices] - right.values[indices])
    })
  nominal <- quantile(draws, c(0.025, 0.975), names = FALSE)
  bonferroni <- quantile(
    draws,
    c(0.05 / (2 * family.size), 1 - 0.05 / (2 * family.size)),
    names = FALSE
    )
  return(tibble(
    difference = mean(left.values) - mean(right.values),
    ci.95.lower = nominal[[1L]], ci.95.upper = nominal[[2L]],
    bonferroni.ci.lower = bonferroni[[1L]],
    bonferroni.ci.upper = bonferroni[[2L]]
    ))
  }


# bootstrap direct population LD differences for each T.C. source and bin
make.ld.population.contrast.tables <- function(
    data, bootstrap.replicates = LD.POPULATION.CONTRAST.BOOTSTRAP.REPLICATES,
    seed = 1L
  ) {
  populations <- unique(c(
    LD.POPULATION.CONTRASTS$simulation.left,
    LD.POPULATION.CONTRASTS$simulation.right
    ))
  validate.ld.paired.ids(
    data %>%
      filter(
        data.type %in% LD.POPULATION.CONTRAST.SOURCES,
        as.character(chrom) %in% SELECTED.CHROMOSOMES,
        pop %in% populations,
        distance_bin_bp %in% LD.PLOT.DISTANCE.BINS
        ),
    c("data.type", "chrom", "distance_bin_bp"), "pop",
    "Paired LD population contrast replicate IDs must match"
    )
  data <- validate.ld.contrast.simulation(
    data, LD.POPULATION.CONTRAST.SOURCES, populations
    )
  set.seed(seed)
  tables <- map_dfr(LD.POPULATION.CONTRAST.SOURCES, function(source) {
    map_dfr(seq_len(nrow(LD.POPULATION.CONTRASTS)), function(index) {
      contrast <- LD.POPULATION.CONTRASTS[index, ]
      map_dfr(LD.PLOT.DISTANCE.BINS, function(distance.bin) {
        values <- data %>%
          filter(
            data.type == source, distance_bin_bp == distance.bin,
            pop %in% c(contrast$simulation.left, contrast$simulation.right)
            )
        left <- values %>% filter(pop == contrast$simulation.left) %>%
          arrange(rep)
        right <- values %>% filter(pop == contrast$simulation.right) %>%
          arrange(rep)
        return(bind_cols(
          tibble(
            data.type = source, distance_bin_bp = distance.bin,
            contrast = contrast$contrast
            ),
          summarize.ld.population.contrast.bootstrap(
            left$mean.r2, right$mean.r2, bootstrap.replicates
            )
          ))
        })
      })
    })
  if (nrow(tables) != length(LD.POPULATION.CONTRAST.SOURCES) *
      nrow(LD.POPULATION.CONTRASTS) * length(LD.PLOT.DISTANCE.BINS)) {
    stop("LD population contrast bootstrap table does not match family size")
    }
  return(tables)
  }


# bootstrap one paired or independent ADX source difference
summarize.ld.simulation.contrast.bootstrap <- function(
    left.values, right.values, bootstrap.replicates, paired,
    family.size = LD.POPULATION.CONTRAST.FAMILY.SIZE
  ) {
  if (length(left.values) != 50L || length(right.values) != 50L ||
      any(!is.finite(left.values)) || any(!is.finite(right.values))) {
    stop("LD simulation contrast bootstraps require 50 finite values")
    }
  draws <- replicate(bootstrap.replicates, {
    if (paired) {
      indices <- sample(seq_along(left.values), length(left.values),
        replace = TRUE)
      mean(left.values[indices] - right.values[indices])
      } else {
      mean(sample(left.values, length(left.values), replace = TRUE)) -
        mean(sample(right.values, length(right.values), replace = TRUE))
      }
    })
  nominal <- quantile(draws, c(0.025, 0.975), names = FALSE)
  bonferroni <- quantile(
    draws,
    c(0.05 / (2 * family.size), 1 - 0.05 / (2 * family.size)),
    names = FALSE
    )
  return(tibble(
    difference = mean(left.values) - mean(right.values),
    ci.95.lower = nominal[[1L]], ci.95.upper = nominal[[2L]],
    bonferroni.ci.lower = bonferroni[[1L]],
    bonferroni.ci.upper = bonferroni[[2L]]
    ))
  }


# bootstrap ADX LD differences among simulation sources for each displayed bin
make.ld.simulation.contrast.tables <- function(
    data, bootstrap.replicates = LD.POPULATION.CONTRAST.BOOTSTRAP.REPLICATES,
    seed = 1L
  ) {
  paired.contrasts <- filter(LD.SIMULATION.CONTRASTS, paired)
  for (index in seq_len(nrow(paired.contrasts))) {
    contrast <- paired.contrasts[index, ]
    validate.ld.paired.ids(
      data %>%
        filter(
          data.type %in% c(contrast$left.source, contrast$right.source),
          as.character(chrom) %in% SELECTED.CHROMOSOMES,
          pop == "ADX", distance_bin_bp %in% LD.PLOT.DISTANCE.BINS
          ),
      c("chrom", "distance_bin_bp"), "data.type",
      "Paired LD simulation contrast replicate IDs must match"
      )
    }
  data <- validate.ld.contrast.simulation(
    data, LD.SIMULATION.CONTRAST.SOURCES, "ADX"
    )
  set.seed(seed)
  tables <- map_dfr(seq_len(nrow(LD.SIMULATION.CONTRASTS)), function(index) {
    contrast <- LD.SIMULATION.CONTRASTS[index, ]
    map_dfr(LD.PLOT.DISTANCE.BINS, function(distance.bin) {
      left <- data %>%
        filter(data.type == contrast$left.source,
          distance_bin_bp == distance.bin) %>%
        arrange(rep)
      right <- data %>%
        filter(data.type == contrast$right.source,
          distance_bin_bp == distance.bin) %>%
        arrange(rep)
      return(bind_cols(
        tibble(distance_bin_bp = distance.bin, contrast = contrast$contrast),
        summarize.ld.simulation.contrast.bootstrap(
          left$mean.r2, right$mean.r2, bootstrap.replicates, contrast$paired
          )
        ))
      })
    })
  if (nrow(tables) != nrow(LD.SIMULATION.CONTRASTS) *
      length(LD.PLOT.DISTANCE.BINS)) {
    stop("LD simulation contrast bootstrap table does not match family size")
    }
  return(tables)
  }


# calculate direct empirical population LD differences for displayed bins
make.ld.population.empirical.contrasts <- function(data) {
  populations <- unique(c(
    LD.POPULATION.CONTRASTS$empirical.left,
    LD.POPULATION.CONTRASTS$empirical.right
    ))
  data <- data %>%
    filter(
      data.type == "Empirical", as.character(chrom) %in% SELECTED.CHROMOSOMES,
      pop %in% populations, distance_bin_bp %in% LD.PLOT.DISTANCE.BINS
      ) %>%
    mutate(
      pop = as.character(pop),
      distance_bin_bp = as.numeric(distance_bin_bp)
      )
  expected <- crossing(
    pop = populations, distance_bin_bp = LD.PLOT.DISTANCE.BINS
    )
  counts <- data %>%
    count(pop, distance_bin_bp, name = "value.count") %>%
    right_join(expected, by = c("pop", "distance_bin_bp")) %>%
    mutate(value.count = replace_na(value.count, 0L)) %>%
    filter(value.count != 1L)
  if (nrow(counts) || any(!is.finite(data$mean.r2))) {
    stop("Empirical LD population contrast inputs require exactly one finite ",
      "value for every selected chromosome, population, and displayed bin")
    }
  return(map_dfr(seq_len(nrow(LD.POPULATION.CONTRASTS)), function(index) {
    contrast <- LD.POPULATION.CONTRASTS[index, ]
    left <- data %>% filter(pop == contrast$empirical.left) %>%
      select(distance_bin_bp, left.value = mean.r2)
    right <- data %>% filter(pop == contrast$empirical.right) %>%
      select(distance_bin_bp, right.value = mean.r2)
    inner_join(left, right, by = "distance_bin_bp") %>%
      transmute(
        distance_bin_bp, contrast = contrast$contrast,
        difference = left.value - right.value
        )
    }))
  }


# select one interval family and significant differences for a contrast plot
prepare.ld.population.contrast.plot.data <- function(
    simulation, interval.type = c("95", "bonferroni"), source
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
  simulation <- simulation %>%
    filter(data.type == source) %>%
    transmute(
      distance_bin_bp, contrast, difference,
      ci.lower = .data[[lower.column]], ci.upper = .data[[upper.column]]
      ) %>%
    mutate(
      contrast = factor(
        contrast, levels = LD.POPULATION.CONTRASTS$contrast
        ),
      point.color = if_else(
        ci.lower > 0 | ci.upper < 0, "red", as.character(contrast)
        )
      )
  return(list(
    simulation = simulation,
    red.markers = filter(simulation, ci.lower > 0 | ci.upper < 0)
    ))
  }


# select one interval family for ADX source-contrast plotting
prepare.ld.simulation.contrast.plot.data <- function(
    simulation, interval.type = c("95", "bonferroni")
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
  simulation <- simulation %>%
    transmute(
      distance_bin_bp, contrast, difference,
      ci.lower = .data[[lower.column]], ci.upper = .data[[upper.column]]
      ) %>%
    mutate(
      contrast = factor(
        contrast, levels = LD.SIMULATION.CONTRASTS$contrast
        ),
      point.color = if_else(
        ci.lower > 0 | ci.upper < 0, "red", as.character(contrast)
        )
      )
  return(list(
    simulation = simulation,
    red.markers = filter(simulation, ci.lower > 0 | ci.upper < 0)
    ))
  }


# construct a simulated LD contrast plot with selected confidence intervals
make.ld.simulation.contrast.plot <- function(
    data, interval.type = c("95", "bonferroni"),
    contrast.type = c("population", "simulation")
  ) {
  interval.type <- match.arg(interval.type)
  contrast.type <- match.arg(contrast.type)
  contrasts <- if (contrast.type == "population") {
    LD.POPULATION.CONTRASTS$contrast
    } else {
    LD.SIMULATION.CONTRASTS$contrast
    }
  colors <- if (contrast.type == "population") {
    PLOT.STYLES$contrast.colors
    } else {
    PLOT.STYLES$simulation.contrast.colors
    }
  active.contrasts <- contrasts[
    contrasts %in% unique(as.character(data$simulation$contrast))
    ]
  title <- paste(
    if (interval.type == "bonferroni") "Bonferroni" else "95%",
    if (contrast.type == "population") {
      "population LD differences"
      } else {
      "ADX simulation LD differences"
      }
    )
  plot <- ggplot(
    data$simulation,
    aes(distance_bin_bp, difference, fill = contrast)
    ) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1) +
    geom_ribbon(aes(ymin = ci.lower, ymax = ci.upper), alpha = 0.2,
      color = NA) +
    geom_line(linewidth = 1) +
    geom_point(
      data = data$simulation,
      aes(
        distance_bin_bp, difference, color = point.color, fill = contrast
        ),
      shape = 21, size = 2.5,
      stroke = LD.CONTRAST.LINEWIDTH, inherit.aes = FALSE
      ) +
    facet_wrap(~contrast, nrow = 1) +
    scale_color_manual(values = c(colors, red = "red"), guide = "none") +
    scale_fill_manual(
      values = colors, breaks = active.contrasts, guide = "none"
      ) +
    scale_x_continuous(
      limits = c(5000, 200000), breaks = LD.PLOT.DISTANCE.BINS
      ) +
    labs(x = "Distance between SNPs (bp)", y = "Difference", title = title) +
    theme()
  return(apply.standard.plot.theme(plot, legend.position = "none"))
  }


# construct a direct empirical LD population-difference plot without intervals
make.ld.population.empirical.contrast.plot <- function(data) {
  data <- data %>%
    mutate(contrast = factor(
      contrast, levels = LD.POPULATION.CONTRASTS$contrast
      ))
  plot <- ggplot(data, aes(distance_bin_bp, difference, color = contrast)) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1) +
    geom_line(linewidth = 1) +
    facet_wrap(~contrast, nrow = 1) +
    scale_color_manual(
      values = PLOT.STYLES$contrast.colors,
      breaks = LD.POPULATION.CONTRASTS$contrast, guide = "none"
      ) +
    scale_x_continuous(
      limits = c(5000, 200000), breaks = LD.PLOT.DISTANCE.BINS
      ) +
    labs(
      x = "Distance between SNPs (bp)", y = "Difference",
      title = "Empirical population LD differences"
      ) +
    theme()
  return(apply.standard.plot.theme(plot, legend.position = "none"))
  }


# build one chromosome-1 bootstrap LD view over the requested distance range
make.bootstrap.ld.plot <- function(data, data.types, view, title = view) {
  source.view <- grepl("all.datatypes.adx.asw", view)
  show.legend <- source.view || view %in% BOOTSTRAP.LEGEND.VIEWS
  plotted <- data %>%
    filter(
      as.character(chrom) == "1",
      as.character(data.type) %in% data.types,
      between(distance_bin_bp, 5000, 200000)
      ) %>%
    filter(
      !source.view |
        (data.type != "Empirical" & pop == "ADX") |
        (data.type == "Empirical" & pop == "ASW")
      )
  active.sources <- order.active.levels(plotted$data.type, SOURCE.LEVELS)
  active.populations <- order.active.levels(
    plotted$pop, POPULATION.LEVELS
    )
  plotted <- plotted %>%
    mutate(
      data.type = factor(as.character(data.type), levels = active.sources),
      pop = factor(as.character(pop), levels = active.populations),
      plot.key = case_when(
        source.view & data.type == "Empirical" ~ "ASW",
        source.view ~ as.character(data.type),
        TRUE ~ as.character(pop)
        ),
      plot.key = factor(
        plot.key,
        levels = if (source.view) {
          c(active.sources, "ASW")[c(active.sources, "ASW") %in% plot.key]
          } else {
          active.populations
          }
        )
      )
  plot <- ggplot(plotted, aes(distance_bin_bp, mean, color = plot.key,
    fill = plot.key, group = interaction(data.type, pop))) +
    scale_x_continuous(
      limits = c(5000, 200000), breaks = LD.PLOT.DISTANCE.BINS
      ) +
    labs(
      x = "Distance between SNPs (bp)", y = expression("Mean " * r^2),
      title = title, color = NULL, fill = NULL
      ) +
    theme()
  if (!grepl("all.lines$", view)) {
    plot <- plot + geom_ribbon(
      data = filter(plotted, data.type != "Empirical"),
      aes(ymin = lower, ymax = upper), alpha = 0.2, color = NA
      )
    }
  source.colors <- c(
    PLOT.STYLES$source.colors,
    ASW = PLOT.STYLES$population.colors[["ASW"]]
    )
  source.labels <- c(PLOT.STYLES$series.labels, ASW = "ASW")
  plot <- plot +
    geom_line(linewidth = 1, show.legend = show.legend) +
    scale_color_manual(
      values = if (source.view) {
        source.colors
        } else {
        PLOT.STYLES$population.colors
        },
      breaks = levels(plotted$plot.key),
      labels = if (source.view) {
        source.labels[levels(plotted$plot.key)]
        } else {
        waiver()
        }
      ) +
    scale_fill_manual(
      values = if (source.view) {
        source.colors
        } else {
        PLOT.STYLES$population.colors
        },
      breaks = levels(plotted$plot.key),
      labels = if (source.view) {
        source.labels[levels(plotted$plot.key)]
        } else {
        waiver()
        }
      ) +
    guides(
      color = if (show.legend) {
        guide_legend(order = 1, nrow = 1, byrow = TRUE)
        } else {
        "none"
        },
      fill = "none"
      )
  if (view == "role.interval") plot <- plot + facet_wrap(~role)
  if (grepl("datatype.interval$", view)) {
    plot <- plot + facet_wrap(
      ~data.type,
      labeller = labeller(data.type = PLOT.STYLES$series.labels)
      )
    }
  return(apply.standard.plot.theme(
    plot, legend.position = if (show.legend) "top" else "none"
    ))
  }


# analysis ----


# read all available chromosomes for all four simulation sources
ld.sim.tc.chromosomes <- read.ld.chromosomes(
  SIM.TC.DATA.DIR, CHROMOSOMES, "Simulation_2T12Consistent"
  )
ld.simDown.tc.chromosomes <- read.ld.chromosomes(
  SIMDOWN.TC.DATA.DIR, CHROMOSOMES,
  "Simulation_2T12Consistent_simDown"
  )
ld.sim.lg.chromosomes <- read.ld.chromosomes(
  SIM.LG.DATA.DIR, CHROMOSOMES, "Simulation_largeGrowth"
  )
ld.simDown.lg.chromosomes <- read.ld.chromosomes(
  SIMDOWN.LG.DATA.DIR, CHROMOSOMES,
  "Simulation_largeGrowth_simDown"
  )
ld.simulation.chromosomes <- bind_rows(
  ld.sim.tc.chromosomes, ld.simDown.tc.chromosomes,
  ld.sim.lg.chromosomes, ld.simDown.lg.chromosomes
  )

# pool every simulation chromosome curve from producer sufficient statistics
ld.simulation.selected <- ld.simulation.chromosomes %>%
  pool.ld.curves(include.chromosome = TRUE)

# read and pool all available empirical chromosomes
ld.empirical.selected <- read.ld.chromosomes(
  EMPIRICAL.DATA.DIR, CHROMOSOMES, "Empirical"
  ) %>%
  pool.ld.curves(include.chromosome = TRUE)
ld.empirical.genome <- read.empirical.ld.genome(EMPIRICAL.DATA.DIR) %>%
  pool.ld.curves(include.chromosome = TRUE)

# summarize curves for the active bootstrap views
ld.summary <- bind_rows(
  ld.simulation.selected, ld.empirical.selected, ld.empirical.genome
  ) %>%
  summarize.ld.curves(SELECTED.CHROMOSOMES)

# bootstrap population and ADX source LD contrasts on displayed bins
ld.population.contrast.tables <- make.ld.population.contrast.tables(
  ld.simulation.selected, LD.POPULATION.CONTRAST.BOOTSTRAP.REPLICATES
  )
ld.simulation.contrast.tables <- make.ld.simulation.contrast.tables(
  ld.simulation.selected, LD.POPULATION.CONTRAST.BOOTSTRAP.REPLICATES
  )
ld.population.contrast.plots <- list(
  tc = list(
    `95` = make.ld.simulation.contrast.plot(
      prepare.ld.population.contrast.plot.data(
        ld.population.contrast.tables, "95", SOURCE.LEVELS[[1L]]
        ),
      "95", "population"
      ),
    bonferroni = make.ld.simulation.contrast.plot(
      prepare.ld.population.contrast.plot.data(
        ld.population.contrast.tables, "bonferroni", SOURCE.LEVELS[[1L]]
        ),
      "bonferroni", "population"
      )
    ),
  tcd = list(
    `95` = make.ld.simulation.contrast.plot(
      prepare.ld.population.contrast.plot.data(
        ld.population.contrast.tables, "95", SOURCE.LEVELS[[2L]]
        ),
      "95", "population"
      ),
    bonferroni = make.ld.simulation.contrast.plot(
      prepare.ld.population.contrast.plot.data(
        ld.population.contrast.tables, "bonferroni", SOURCE.LEVELS[[2L]]
        ),
      "bonferroni", "population"
      )
    ),
  empirical = make.ld.population.empirical.contrast.plot(
    make.ld.population.empirical.contrasts(ld.empirical.selected)
    )
  )
ld.simulation.contrast.plots <- list(
  `95` = make.ld.simulation.contrast.plot(
    prepare.ld.simulation.contrast.plot.data(
      ld.simulation.contrast.tables, "95"
      ),
    "95", "simulation"
    ),
  bonferroni = make.ld.simulation.contrast.plot(
    prepare.ld.simulation.contrast.plot.data(
      ld.simulation.contrast.tables, "bonferroni"
      ),
    "bonferroni", "simulation"
    )
  )

# save TCD/1kG and all-datatype ADX/ASW chromosome-1 bootstrap LD views
ld.bootstrap.tcd.1kg.all.lines <- make.bootstrap.ld.plot(
  ld.summary, c("Simulation_2T12Consistent_simDown", "Empirical"),
  "all.lines", "LD decay: chromosome 1 TCD and 1kG"
  )
ld.bootstrap.tcd.1kg.role.interval <- make.bootstrap.ld.plot(
  ld.summary, c("Simulation_2T12Consistent_simDown", "Empirical"),
  "role.interval", "LD decay: chromosome 1 TCD and 1kG by role"
  )
ld.bootstrap.tcd.1kg.datatype.interval <- make.bootstrap.ld.plot(
  ld.summary, c("Simulation_2T12Consistent_simDown", "Empirical"),
  "datatype.interval", "LD decay: chromosome 1 TCD and 1kG by source"
  )
ld.bootstrap.all.datatypes.adx.asw.all.lines <- make.bootstrap.ld.plot(
  ld.summary, SOURCE.LEVELS, "all.datatypes.adx.asw.all.lines",
  "LD decay: chromosome 1 all ADX sources and ASW"
  )
ld.bootstrap.all.datatypes.adx.asw.datatype.interval <- make.bootstrap.ld.plot(
  ld.summary, SOURCE.LEVELS, "all.datatypes.adx.asw.datatype.interval",
  "LD decay: chromosome 1 all ADX sources and ASW by source"
  )
dir.create(OUTPUT.DIR, recursive = TRUE, showWarnings = FALSE)
saveRDS(ld.bootstrap.tcd.1kg.all.lines, file.path(
  OUTPUT.DIR, "ld.bootstrap.tcd.1kg.all.lines.rds"
  ))
saveRDS(ld.bootstrap.tcd.1kg.role.interval, file.path(
  OUTPUT.DIR, "ld.bootstrap.tcd.1kg.role.interval.rds"
  ))
saveRDS(ld.bootstrap.tcd.1kg.datatype.interval, file.path(
  OUTPUT.DIR, "ld.bootstrap.tcd.1kg.datatype.interval.rds"
  ))
saveRDS(ld.bootstrap.all.datatypes.adx.asw.all.lines, file.path(
  OUTPUT.DIR, "ld.bootstrap.all.datatypes.adx.asw.all.lines.rds"
  ))
saveRDS(ld.bootstrap.all.datatypes.adx.asw.datatype.interval, file.path(
  OUTPUT.DIR, "ld.bootstrap.all.datatypes.adx.asw.datatype.interval.rds"
  ))
saveRDS(ld.population.contrast.plots$tc$`95`, file.path(
  OUTPUT.DIR, "ld.population.tc.95.rds"
  ))
saveRDS(ld.population.contrast.plots$tc$bonferroni, file.path(
  OUTPUT.DIR, "ld.population.tc.bonferroni.rds"
  ))
saveRDS(ld.population.contrast.plots$tcd$`95`, file.path(
  OUTPUT.DIR, "ld.population.tcd.95.rds"
  ))
saveRDS(ld.population.contrast.plots$tcd$bonferroni, file.path(
  OUTPUT.DIR, "ld.population.tcd.bonferroni.rds"
  ))
saveRDS(ld.population.contrast.plots$empirical, file.path(
  OUTPUT.DIR, "ld.population.empirical.rds"
  ))
saveRDS(ld.simulation.contrast.plots$`95`, file.path(
  OUTPUT.DIR, "ld.simulation.contrast.95.rds"
  ))
saveRDS(ld.simulation.contrast.plots$bonferroni, file.path(
  OUTPUT.DIR, "ld.simulation.contrast.bonferroni.rds"
  ))

print(ld.bootstrap.tcd.1kg.all.lines)
# print(ld.bootstrap.tcd.1kg.role.interval)
print(ld.bootstrap.tcd.1kg.datatype.interval)
print(ld.bootstrap.all.datatypes.adx.asw.all.lines)
# print(ld.bootstrap.all.datatypes.adx.asw.datatype.interval)
# print(ld.population.contrast.plots$tc$`95`)
print(ld.population.contrast.plots$tc$bonferroni)
# print(ld.population.contrast.plots$tcd$`95`)
print(ld.population.contrast.plots$tcd$bonferroni)
print(ld.population.contrast.plots$empirical)
# print(ld.simulation.contrast.plots$`95`)
print(ld.simulation.contrast.plots$bonferroni)

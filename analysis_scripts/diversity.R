# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology
# Mooney Lab
# ___
# diversity.R
# ______________________________________________________________________________

# pattern: Mixed (unavoidable)
# Reason: This analysis script combines Parquet I/O, bootstrap calculations,
# and plot persistence.


# set up ----
library(tidyverse)
library(nanoparquet)


SIM.TC.DATA.DIR <- "~/scratch/OOA_NAAdmixture_2T12Consistent/stats"
SIMDOWN.TC.DATA.DIR <- "~/scratch/OOA_NAAdmixture_2T12ConsistentOnekgDownsample/stats"
SIM.LG.DATA.DIR <- "~/scratch/OOA_NAAdmixture_largeGrowth/stats"
SIMDOWN.LG.DATA.DIR <- "~/scratch/OOA_NAAdmixture_largeGrowthOnekgDownsample/stats"
EMPIRICAL.DATA.DIR <- "~/scratch/OOA_NAAdmixture_1kG/stats"
OUTPUT.DIR <- "/home1/karatas/proj/OOA_NAAdmixture_data"
CHROMOSOMES <- as.character(1:22)
SELECTED.CHROMOSOMES <- c("1", "10", "20")
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
PLOT.BASE.SIZE <- 24
CATEGORICAL.BAR.DODGE <- 0.9
CATEGORICAL.BAR.WIDTH <- 0.8
CATEGORICAL.BAR.LINEWIDTH <- 1
DIVERSITY.POPULATION.CONTRAST.FAMILY.SIZE <- 66L
DIVERSITY.POPULATION.CONTRAST.BOOTSTRAP.REPLICATES <- 1000L
DIVERSITY.POPULATION.CONTRAST.SOURCES <- SOURCE.LEVELS[seq_len(4L)]
DIVERSITY.POPULATION.CONTRAST.SOURCE.GROUPS <- list(
  tc.tcd = SOURCE.LEVELS[c(1L, 2L)],
  lg.lgd = SOURCE.LEVELS[c(3L, 4L)]
  )
DIVERSITY.POPULATION.CONTRASTS <- tribble(
  ~contrast, ~simulation.left, ~simulation.right,
  ~empirical.left, ~empirical.right,
  "AFR-ADX", "AFR", "ADX", "YRI", "ASW",
  "AFR-EUR", "AFR", "EUR", "YRI", "CEU",
  "ADX-EUR", "ADX", "EUR", "ASW", "CEU"
  )
DIVERSITY.SIMULATION.CONTRASTS <- tribble(
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
  fill.colors = c(
    "TC AFR" = "#9BD5F2", "TC ADX" = "#9A83CE",
    "TC EUR" = "#FBB4AE",
    "TC D. AFR" = "#56B4E9",
    "TC D. ADX" = "#6F55B5",
    "TC D. EUR" = "#FB8072",
    "LG ADX" = "#32146F",
    "LG D. ADX" = "#4B1FA8",
    "empirical YRI" = "#EEC4DC", "empirical ASW" = "#E44B8D",
    "empirical CEU" = "#BB437E"
    ),
  fill.labels = c(
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
    ),
  empirical.colors = c(
    YRI = "#EEC4DC", ASW = "#E44B8D", CEU = "#BB437E"
    ),
  contrast.colors = c(
    `AFR-ADX` = "#BFFBFF", `AFR-EUR` = "#16ACBD",
    `ADX-EUR` = "#00606F"
    ),
  contrast.labels = c(
    `AFR-ADX` = "AFR - ADX", `AFR-EUR` = "AFR - EUR",
    `ADX-EUR` = "ADX - EUR"
    ),
  contrast.linetypes = c(
    `AFR-ADX` = "solid", `AFR-EUR` = "longdash",
    `ADX-EUR` = "dotted"
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


# summarize complete simulation replicates with direct percentile intervals
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


# retain source-specific populations before any replicate summaries
apply.diversity.source.contract <- function(data) {
  retained <- data %>%
    filter(
      (data.type %in% c(
        "Simulation_2T12Consistent", "Simulation_2T12Consistent_simDown"
        ) & pop %in% c("AFR", "ADX", "EUR")) |
        (data.type %in% c(
          "Simulation_largeGrowth", "Simulation_largeGrowth_simDown"
        ) & pop %in% c("AFR", "ADX", "EUR")) |
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
    stop("Diversity data contain an unsupported population label")
    }

  return(data)
  }


# validate complete intergenic simulation values for population contrasts
validate.diversity.population.contrast.simulation <- function(
    data, chromosomes = CHROMOSOMES
  ) {
  required.columns <- c(
    "data.type", "rep", "chrom", "pop", "stat", "value", "mask"
    )
  missing.columns <- setdiff(required.columns, names(data))
  if (length(missing.columns)) {
    stop("Population contrast data are missing columns: ",
      paste(missing.columns, collapse = ", "))
    }
  chromosomes <- as.character(chromosomes)
  sources <- DIVERSITY.POPULATION.CONTRAST.SOURCES
  populations <- unique(c(
    DIVERSITY.POPULATION.CONTRASTS$simulation.left,
    DIVERSITY.POPULATION.CONTRASTS$simulation.right
    ))
  data <- data %>%
    filter(
      data.type %in% sources, chrom %in% chromosomes,
      pop %in% populations, stat %in% c("pi", "theta"),
      mask == "Intergenic"
      ) %>%
    mutate(
      data.type = as.character(data.type), chrom = as.character(chrom),
      pop = as.character(pop), stat = as.character(stat)
      )
  if (any(!is.finite(data$value))) {
    stop("Population contrast simulation values must be finite")
    }
  group.counts <- data %>%
    count(data.type, chrom, stat, pop, name = "row.count") %>%
    filter(row.count != 50L)
  if (nrow(group.counts)) {
    stop("Population contrast simulation inputs require exactly one finite ",
      "value for every source, replicate, chromosome, population, and ",
      "statistic")
    }
  paired.ids <- data %>%
    group_by(data.type, chrom, stat, pop) %>%
    summarise(rep.ids = list(sort(unique(rep))), .groups = "drop") %>%
    pivot_wider(names_from = pop, values_from = rep.ids)
  for (index in seq_len(nrow(paired.ids))) {
    values <- unlist(paired.ids[index, populations], recursive = FALSE)
    if (length(values) != length(populations) ||
        !all(vapply(values, identical, logical(1), values[[1L]]))) {
      stop("Population contrast simulation replicate IDs must match")
      }
    }
  expected <- crossing(
    data.type = sources, rep = seq_len(50L), chrom = chromosomes,
    pop = populations, stat = c("pi", "theta")
    )
  counts <- expected %>%
    left_join(
      data %>% count(data.type, rep, chrom, pop, stat, name = "value.count"),
      by = c("data.type", "rep", "chrom", "pop", "stat")
      ) %>%
    mutate(value.count = replace_na(value.count, 0L)) %>%
    filter(value.count != 1L)
  if (nrow(counts)) {
    stop("Population contrast simulation inputs require exactly one finite ",
      "value for every source, replicate, chromosome, population, and ",
      "statistic")
    }
  unexpected <- anti_join(
    data,
    expected,
    by = c("data.type", "rep", "chrom", "pop", "stat")
    )
  if (nrow(unexpected)) {
    stop("Population contrast simulation replicate IDs must match")
    }
  return(data)
  }


# summarize direct bootstrap intervals for one population difference
summarize.diversity.population.contrast.bootstrap <- function(
    left.values, right.values, bootstrap.replicates,
    family.size = DIVERSITY.POPULATION.CONTRAST.FAMILY.SIZE
  ) {
  if (length(left.values) != 50L || length(right.values) != 50L ||
      any(!is.finite(left.values)) || any(!is.finite(right.values))) {
    stop("Population contrast bootstraps require 50 finite paired values")
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


# bootstrap paired population differences for all sources and chromosomes
make.diversity.population.contrast.tables <- function(
    data, chromosomes = CHROMOSOMES,
    bootstrap.replicates = DIVERSITY.POPULATION.CONTRAST.BOOTSTRAP.REPLICATES,
    seed = 1L
  ) {
  chromosomes <- as.character(chromosomes)
  data <- validate.diversity.population.contrast.simulation(data, chromosomes)
  if (!is.numeric(bootstrap.replicates) || length(bootstrap.replicates) != 1L ||
      !is.finite(bootstrap.replicates) || bootstrap.replicates < 1L ||
      bootstrap.replicates %% 1L != 0) {
    stop(
      "Population contrast bootstrap replicate count must be a positive ",
      "integer"
      )
    }
  set.seed(seed)
  tables <- map_dfr(DIVERSITY.POPULATION.CONTRAST.SOURCES, function(source) {
    map_dfr(c("pi", "theta"), function(statistic) {
      map_dfr(chromosomes, function(chromosome) {
        map_dfr(seq_len(nrow(DIVERSITY.POPULATION.CONTRASTS)), function(index) {
          contrast <- DIVERSITY.POPULATION.CONTRASTS[index, ]
          values <- data %>%
            filter(
              data.type == source, stat == statistic, chrom == chromosome,
              pop %in% c(contrast$simulation.left, contrast$simulation.right)
              ) %>%
            arrange(pop, rep)
          left.values <- values %>%
            filter(pop == contrast$simulation.left) %>% arrange(rep) %>%
            pull(value)
          right.values <- values %>%
            filter(pop == contrast$simulation.right) %>% arrange(rep) %>%
            pull(value)
          return(bind_cols(
            tibble(
              data.type = source, stat = statistic, chrom = chromosome,
              contrast = contrast$contrast
              ),
            summarize.diversity.population.contrast.bootstrap(
              left.values, right.values, bootstrap.replicates
              )
            ))
          })
        })
      })
    })
  expected.rows <- length(DIVERSITY.POPULATION.CONTRAST.SOURCES) * 2L *
    DIVERSITY.POPULATION.CONTRAST.FAMILY.SIZE
  if (nrow(tables) != expected.rows) {
    stop("Population contrast bootstrap table does not match its family size")
    }
  return(tables)
  }


# validate complete ADX simulation values for source contrasts
validate.diversity.simulation.contrast.data <- function(
    data, chromosomes = CHROMOSOMES
  ) {
  required.columns <- c(
    "data.type", "rep", "chrom", "pop", "stat", "value", "mask"
    )
  missing.columns <- setdiff(required.columns, names(data))
  if (length(missing.columns)) {
    stop("Simulation contrast data are missing columns: ",
      paste(missing.columns, collapse = ", "))
    }
  chromosomes <- as.character(chromosomes)
  data <- data %>%
    filter(
      data.type %in% DIVERSITY.POPULATION.CONTRAST.SOURCES,
      chrom %in% chromosomes, pop == "ADX", stat %in% c("pi", "theta"),
      mask == "Intergenic"
      ) %>%
    mutate(data.type = as.character(data.type), chrom = as.character(chrom))
  if (any(!is.finite(data$value))) {
    stop("Simulation contrast values must be finite")
    }
  paired.contrasts <- filter(DIVERSITY.SIMULATION.CONTRASTS, paired)
  for (index in seq_len(nrow(paired.contrasts))) {
    contrast <- paired.contrasts[index, ]
    ids <- data %>%
      filter(data.type %in% c(contrast$left.source, contrast$right.source)) %>%
      group_by(chrom, stat, data.type) %>%
      summarise(rep.ids = list(sort(rep)), .groups = "drop") %>%
      pivot_wider(names_from = data.type, values_from = rep.ids)
    for (row in seq_len(nrow(ids))) {
      if (!identical(ids[[contrast$left.source]][[row]],
          ids[[contrast$right.source]][[row]])) {
        stop("Paired simulation contrast replicate IDs must match")
        }
      }
    }
  expected <- crossing(
    data.type = DIVERSITY.POPULATION.CONTRAST.SOURCES, rep = seq_len(50L),
    chrom = chromosomes, pop = "ADX", stat = c("pi", "theta")
    )
  counts <- data %>%
    count(data.type, chrom, stat, rep, name = "value.count") %>%
    filter(value.count != 1L)
  replicate.counts <- data %>%
    count(data.type, chrom, stat, name = "replicate.count") %>%
    filter(replicate.count != 50L)
  expected.counts <- expected %>%
    left_join(
      data %>% count(data.type, rep, chrom, pop, stat, name = "value.count"),
      by = c("data.type", "rep", "chrom", "pop", "stat")
      ) %>%
    mutate(value.count = replace_na(value.count, 0L)) %>%
    filter(value.count != 1L)
  if (nrow(counts) || nrow(replicate.counts) || nrow(expected.counts)) {
    stop("Simulation contrast inputs require exactly one finite value for ",
      "every source, replicate, chromosome, and statistic")
    }
  return(data)
  }


# bootstrap one paired or independent ADX source difference
summarize.diversity.simulation.contrast.bootstrap <- function(
    left.values, right.values, bootstrap.replicates, paired,
    family.size = DIVERSITY.POPULATION.CONTRAST.FAMILY.SIZE
  ) {
  if (length(left.values) != 50L || length(right.values) != 50L ||
      any(!is.finite(left.values)) || any(!is.finite(right.values))) {
    stop("Simulation contrast bootstraps require 50 finite values")
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


# bootstrap ADX source contrasts with their configured replicate relationship
make.diversity.simulation.contrast.tables <- function(
    data, chromosomes = CHROMOSOMES,
    bootstrap.replicates = DIVERSITY.POPULATION.CONTRAST.BOOTSTRAP.REPLICATES,
    seed = 1L
  ) {
  chromosomes <- as.character(chromosomes)
  data <- validate.diversity.simulation.contrast.data(data, chromosomes)
  if (!is.numeric(bootstrap.replicates) || length(bootstrap.replicates) != 1L ||
      !is.finite(bootstrap.replicates) || bootstrap.replicates < 1L ||
      bootstrap.replicates %% 1L != 0) {
    stop("Simulation contrast bootstrap replicate count must be a positive ",
      "integer")
    }
  set.seed(seed)
  tables <- map_dfr(c("pi", "theta"), function(statistic) {
    map_dfr(chromosomes, function(chromosome) {
      map_dfr(seq_len(nrow(DIVERSITY.SIMULATION.CONTRASTS)), function(index) {
        contrast <- DIVERSITY.SIMULATION.CONTRASTS[index, ]
        left <- data %>%
          filter(
            data.type == contrast$left.source, stat == statistic,
            chrom == chromosome
            ) %>%
          arrange(rep)
        right <- data %>%
          filter(
            data.type == contrast$right.source, stat == statistic,
            chrom == chromosome
            ) %>%
          arrange(rep)
        if (contrast$paired && !identical(left$rep, right$rep)) {
          stop("Paired simulation contrast replicate IDs must match")
          }
        return(bind_cols(
          tibble(stat = statistic, chrom = chromosome,
            contrast = contrast$contrast),
          summarize.diversity.simulation.contrast.bootstrap(
            left$value, right$value, bootstrap.replicates, contrast$paired
            )
          ))
        })
      })
    })
  expected.rows <- nrow(DIVERSITY.SIMULATION.CONTRASTS) * 2L *
    length(chromosomes)
  if (nrow(tables) != expected.rows) {
    stop("Simulation contrast bootstrap table does not match its family size")
    }
  return(tables)
  }


# calculate direct empirical chromosome and genome population differences
make.diversity.population.empirical.contrasts <- function(
    data, chromosomes = CHROMOSOMES
  ) {
  chromosomes <- as.character(chromosomes)
  data <- data %>%
    filter(
      data.type == "Empirical", chrom %in% c(chromosomes, "all"),
      pop %in% c("YRI", "ASW", "CEU"), stat %in% c("pi", "theta")
      ) %>%
    mutate(chrom = as.character(chrom), pop = as.character(pop))
  expected <- crossing(
    chrom = c(chromosomes, "all"), pop = c("YRI", "ASW", "CEU"),
    stat = c("pi", "theta")
    )
  counts <- expected %>%
    left_join(data %>% count(chrom, pop, stat, name = "value.count"),
      by = c("chrom", "pop", "stat")) %>%
    mutate(value.count = replace_na(value.count, 0L)) %>%
    filter(value.count != 1L)
  if (nrow(counts) || any(!is.finite(data$value))) {
    stop(
      "Empirical population contrast inputs require exactly one finite value"
      )
    }
  return(map_dfr(
    seq_len(nrow(DIVERSITY.POPULATION.CONTRASTS)), function(index) {
      contrast <- DIVERSITY.POPULATION.CONTRASTS[index, ]
      left <- data %>% filter(pop == contrast$empirical.left) %>%
        select(chrom, stat, left.value = value)
      right <- data %>% filter(pop == contrast$empirical.right) %>%
        select(chrom, stat, right.value = value)
      inner_join(left, right, by = c("chrom", "stat")) %>%
        transmute(chrom, stat, contrast = contrast$contrast,
          difference = left.value - right.value)
      }
    ))
  }


# prepare simulated source-specific values for one contrast interval type
prepare.diversity.population.contrast.plot.data <- function(
    simulation, interval.type = c("95", "bonferroni"),
    chromosomes = CHROMOSOMES,
    sources = DIVERSITY.POPULATION.CONTRAST.SOURCES
  ) {
  interval.type <- match.arg(interval.type)
  chromosomes <- as.character(chromosomes)
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
  contrasts <- DIVERSITY.POPULATION.CONTRASTS$contrast
  simulation <- simulation %>%
    filter(
      data.type %in% sources,
      chrom %in% chromosomes,
      contrast %in% contrasts
      ) %>%
    transmute(
      data.type, stat, chrom, contrast, difference,
      ci.lower = .data[[lower.column]], ci.upper = .data[[upper.column]]
      )
  simulation <- simulation %>%
    mutate(
      data.type = factor(
        data.type, levels = sources
        ),
      chrom = factor(chrom, levels = chromosomes),
      contrast = factor(contrast, levels = contrasts),
      stat = factor(stat, levels = c("pi", "theta")),
      plot.x = as.numeric(chrom) + (as.numeric(contrast) - 2) * 0.25
      )
  markers <- simulation %>%
    group_by(data.type, stat) %>%
    mutate(
      simulated.outside = ci.lower > 0 | ci.upper < 0,
      point.color = if_else(
        simulated.outside, "red", "black"
        )
      ) %>%
    ungroup()
  return(list(
    simulation = markers,
    red.markers = filter(markers, simulated.outside)
    ))
  }


# prepare empirical chromosome points and genome reference lines separately
prepare.diversity.population.empirical.contrast.plot.data <- function(
    empirical, chromosomes = CHROMOSOMES
  ) {
  chromosomes <- as.character(chromosomes)
  contrasts <- DIVERSITY.POPULATION.CONTRASTS$contrast
  chromosome <- empirical %>%
    filter(chrom %in% chromosomes, contrast %in% contrasts) %>%
    mutate(
      chrom = factor(chrom, levels = chromosomes),
      contrast = factor(contrast, levels = contrasts),
      stat = factor(stat, levels = c("pi", "theta")),
      plot.x = as.numeric(chrom) + (as.numeric(contrast) - 2) * 0.25
      )
  genome <- empirical %>%
    filter(chrom == "all", contrast %in% contrasts) %>%
    mutate(
      contrast = factor(contrast, levels = contrasts),
      stat = factor(stat, levels = c("pi", "theta"))
      )
  return(list(chromosome = chromosome, genome = genome))
  }


# prepare one ADX simulation-contrast interval type for plotting
prepare.diversity.simulation.contrast.plot.data <- function(
    simulation, interval.type = c("95", "bonferroni"),
    chromosomes = CHROMOSOMES
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
  contrasts <- DIVERSITY.SIMULATION.CONTRASTS$contrast
  simulation <- simulation %>%
    filter(chrom %in% as.character(chromosomes), contrast %in% contrasts) %>%
    transmute(
      stat, chrom, contrast, difference,
      ci.lower = .data[[lower.column]], ci.upper = .data[[upper.column]]
      ) %>%
    mutate(
      chrom = factor(chrom, levels = as.character(chromosomes)),
      contrast = factor(contrast, levels = contrasts),
      stat = factor(stat, levels = c("pi", "theta")),
      plot.x = as.numeric(chrom) + (as.numeric(contrast) - 2) * 0.25
      ) %>%
    group_by(stat) %>%
    mutate(
      simulated.outside = ci.lower > 0 | ci.upper < 0,
      point.color = if_else(
        simulated.outside, "red", as.character(contrast)
        )
      ) %>%
    ungroup()
  return(list(
    simulation = simulation,
    red.markers = filter(simulation, simulated.outside)
    ))
  }


# construct one simulated diversity contrast plot
make.diversity.contrast.plot <- function(
    data, interval.type = c("95", "bonferroni"),
    contrast.type = c("population", "simulation")
  ) {
  interval.type <- match.arg(interval.type)
  contrast.type <- match.arg(contrast.type)
  simulation <- data$simulation
  contrasts <- if (contrast.type == "population") {
    DIVERSITY.POPULATION.CONTRASTS$contrast
    } else {
    DIVERSITY.SIMULATION.CONTRASTS$contrast
    }
  contrast.colors <- if (contrast.type == "population") {
    PLOT.STYLES$contrast.colors
    } else {
    PLOT.STYLES$simulation.contrast.colors
    }
  active.contrasts <- contrasts[
    contrasts %in% unique(as.character(simulation$contrast))
    ]
  plot <- ggplot(simulation, aes(plot.x, difference)) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1)
  facet.layer <- if (contrast.type == "population") {
    facet_grid(
      rows = vars(data.type, stat), scales = "free_y",
      labeller = labeller(
        data.type = SOURCE.LABELS, stat = c(pi = "π", theta = "θ[w]")
        )
      )
    } else {
    facet_grid(
      rows = vars(stat), scales = "free_y",
      labeller = labeller(stat = c(pi = "π", theta = "θ[w]"))
      )
    }
  plot <- plot +
    geom_errorbar(
      aes(ymin = ci.lower, ymax = ci.upper), width = 0, color = "black",
      linewidth = CATEGORICAL.BAR.LINEWIDTH, show.legend = FALSE
      ) +
    geom_col(
      aes(fill = contrast, color = point.color),
      position = position_dodge(width = CATEGORICAL.BAR.DODGE),
      width = CATEGORICAL.BAR.WIDTH,
      linewidth = CATEGORICAL.BAR.LINEWIDTH
      ) +
    facet.layer +
    scale_color_manual(values = c(black = "black", red = "red"),
      guide = "none") +
    scale_fill_manual(
      values = contrast.colors, breaks = active.contrasts,
      labels = if (contrast.type == "population") {
        PLOT.STYLES$contrast.labels[active.contrasts]
        } else {
        active.contrasts
        }
      ) +
    scale_linetype_manual(
      values = PLOT.STYLES$contrast.linetypes, guide = "none"
      ) +
    guides(
      fill = guide_legend(
        override.aes = list(
          shape = 22, fill = unname(contrast.colors[active.contrasts]),
          color = "black"
          )
        )
      ) +
    scale_x_continuous(
      breaks = seq_along(CHROMOSOMES), labels = CHROMOSOMES
      ) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.17))) +
    labs(
      x = "Chromosome", y = "Difference", color = NULL,
      title = paste(
        if (interval.type == "bonferroni") "Bonferroni" else "95%",
        if (contrast.type == "population") {
          "population diversity differences"
          } else {
          "ADX simulation diversity differences"
          }
        )
      ) +
    theme(legend.title = element_blank())
  return(apply.standard.plot.theme(plot))
  }


# construct an empirical-only diversity population contrast plot
make.diversity.population.empirical.contrast.plot <- function(data) {
  contrasts <- DIVERSITY.POPULATION.CONTRASTS$contrast
  active.contrasts <- contrasts[
    contrasts %in% unique(as.character(data$chromosome$contrast))
    ]
  plot <- ggplot(data$chromosome, aes(plot.x, difference)) +
    geom_hline(
      data = data$genome,
      aes(yintercept = difference, color = contrast, linetype = contrast),
      alpha = 0.5, linewidth = CATEGORICAL.BAR.LINEWIDTH,
      show.legend = FALSE
      ) +
    geom_col(
      aes(plot.x, difference, fill = contrast),
      position = position_dodge(width = CATEGORICAL.BAR.DODGE),
      width = CATEGORICAL.BAR.WIDTH, color = "black",
      linewidth = CATEGORICAL.BAR.LINEWIDTH, inherit.aes = FALSE
      ) +
    facet_grid(
      rows = vars(stat), scales = "free_y",
      labeller = labeller(stat = c(pi = "π", theta = "θ[w]"))
      ) +
    scale_color_manual(
      values = PLOT.STYLES$contrast.colors, breaks = active.contrasts,
      labels = PLOT.STYLES$contrast.labels[active.contrasts]
      ) +
    scale_fill_manual(values = PLOT.STYLES$contrast.colors, guide = "none") +
    scale_linetype_manual(
      values = PLOT.STYLES$contrast.linetypes, guide = "none"
      ) +
    guides(
      color = guide_legend(
        override.aes = list(
          shape = 22,
          fill = unname(PLOT.STYLES$contrast.colors[active.contrasts]),
          color = "black"
          )
        )
      ) +
    scale_x_continuous(
      breaks = seq_along(CHROMOSOMES), labels = CHROMOSOMES
      ) +
    labs(
      x = "Chromosome", y = "Difference", color = NULL,
      title = "Empirical population diversity differences"
      ) +
    theme(legend.title = element_blank())
  return(apply.standard.plot.theme(plot))
  }


# standardize one diversity table to the shared analysis schema
normalize.diversity.table <- function(
    data, data.type.input, mask.input = NA_character_
  ) {
  if (!"span" %in% names(data)) {
    stop("Diversity data are missing the required span column")
    }
  data <- data %>%
    mutate(
      rep = as.numeric(rep),
      chrom = as.character(chrom),
      pop = as.character(pop),
      stat = as.character(stat),
      value = as.numeric(value),
      span = as.numeric(span),
      data.type = data.type.input,
      mask = mask.input
      ) %>%
    apply.diversity.source.contract() %>%
    add.population.roles()
  if (any(!is.finite(data$span) | data$span <= 0)) {
    stop("Diversity spans must be positive and finite")
    }

  return(data)
  }


# read chromosome-labelled diversity files for one source
read.diversity.chromosomes <- function(
    data.directory, file.family, chromosomes, data.type.input,
    mask.input = NA_character_
  ) {
  paths <- file.path(
    path.expand(data.directory),
    str_replace(file.family, fixed("{chrom}"), chromosomes)
    )
  missing <- !file.exists(paths)
  if (any(missing)) {
    warning(
      paste0(
        data.type.input, " ", file.family,
        " is unavailable for chromosomes: ",
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
    return(normalize.diversity.table(
      table, data.type.input, mask.input
      ))
    })

  return(data)
  }


# read one genome-wide empirical diversity file
read.diversity.genome <- function(
    data.directory, file.name, mask.input
  ) {
  data <- read_parquet(file.path(
    path.expand(data.directory), file.name
    ))
  data$chrom <- "all"
  data <- normalize.diversity.table(
    data, "Empirical", mask.input
    )

  return(data)
  }


# calculate replicate means and standard deviations for simulations
summarize.simulation.diversity <- function(data) {
  if (!"role" %in% names(data)) data <- add.population.roles(data)
  replicate.summary <- data %>%
    filter(stat %in% c("pi", "theta")) %>%
    group_by(data.type, rep, pop, role, stat, chrom, mask, span) %>%
    summarise(
      value = mean(value),
      .groups = "drop"
      )
  summary <- summarize.bootstrap.interval(
    replicate.summary,
    c(
      "data.type", "rep", "pop", "role", "stat", "chrom", "mask",
      "span"
      ),
    "value"
    )
  return(summary)
  }


# duplicate simulation summaries across empirical mask comparisons
duplicate.simulation.masks <- function(data) {
  original <- data %>%
    filter(!str_detect(as.character(data.type), "simDown$")) %>%
    select(-mask) %>%
    crossing(mask = c("Intergenic", "Full callable"))
  simDown <- data %>%
    filter(str_detect(as.character(data.type), "simDown$"))
  duplicated <- bind_rows(original, simDown)

  return(duplicated)
  }


# assemble selected-chromosome points and empirical genome references
build.diversity.plot.data <- function(
    simulation.summary, empirical.chromosome, empirical.genome,
    chromosomes
  ) {
  if (!"role" %in% names(empirical.chromosome)) {
    empirical.chromosome <- add.population.roles(empirical.chromosome)
    }
  if (!"role" %in% names(empirical.genome)) {
    empirical.genome <- add.population.roles(empirical.genome)
    }
  if (!"data.type" %in% names(empirical.chromosome)) {
    empirical.chromosome$data.type <- "Empirical"
    }
  simulation.points <- simulation.summary %>%
    filter(chrom %in% chromosomes) %>%
    duplicate.simulation.masks() %>%
    transmute(
      data.type, pop, role, stat, chrom, mask, span,
      estimate = mean, lower, upper, replicate.count
      )
  empirical.points <- empirical.chromosome %>%
    filter(chrom %in% chromosomes, stat %in% c("pi", "theta")) %>%
    transmute(
      data.type, pop, role, stat, chrom, mask, span,
      estimate = value, lower = NA_real_, upper = NA_real_, replicate.count = 1L
      )
  points <- bind_rows(simulation.points, empirical.points) %>%
    mutate(
      chrom = factor(chrom, levels = chromosomes),
      pop = factor(pop, levels = POPULATION.LEVELS),
      data.type = factor(
        data.type, levels = SOURCE.LEVELS
        ),
      mask = factor(mask, levels = c("Intergenic", "Full callable")),
      stat = factor(stat, levels = c("pi", "theta")),
      fill.key = factor(case_when(
        data.type == "Simulation_2T12Consistent" ~ paste("TC", pop),
        data.type == "Simulation_largeGrowth" ~ "LG ADX",
        data.type == "Simulation_2T12Consistent_simDown" ~
          paste("TC D.", pop),
        data.type == "Simulation_largeGrowth_simDown" ~ "LG D. ADX",
        data.type == "Empirical" ~ paste("empirical", pop)
        ), levels = names(PLOT.STYLES$fill.colors))
      )
  genome.lines <- empirical.genome %>%
    filter(chrom == "all", stat %in% c("pi", "theta")) %>%
    transmute(pop, role, stat, mask, span, estimate = value) %>%
    mutate(
      pop = factor(pop, levels = POPULATION.LEVELS),
      mask = factor(mask, levels = c("Intergenic", "Full callable")),
      stat = factor(stat, levels = c("pi", "theta"))
      )

  return(list(points = points, genome.lines = genome.lines))
  }


# retain one configured diversity view and its eligible empirical references
filter.diversity.plot.view <- function(
    points, genome.lines, data.types, tag
  ) {
  if (!tag %in% names(PLOT.CONFIGS)) {
    stop("Unsupported diversity plot tag: ", tag)
    }
  if (!identical(data.types, PLOT.CONFIGS[[tag]])) {
    stop("Diversity data types do not match the configured tag")
    }
  view.points <- points %>%
    filter(as.character(data.type) %in% data.types) %>%
    filter(
      tag != "all.datatypes.adx.asw" |
        (data.type != "Empirical" & pop == "ADX") |
        (data.type == "Empirical" & pop == "ASW")
      )
  active.sources <- order.active.levels(
    view.points$data.type, SOURCE.LEVELS
    )
  active.populations <- order.active.levels(
    view.points$pop, POPULATION.LEVELS
    )
  view.points <- view.points %>%
    mutate(
      data.type = factor(as.character(data.type), levels = active.sources),
      pop = factor(as.character(pop), levels = active.populations),
      fill.key = factor(
        as.character(fill.key),
        levels = order.active.levels(
          fill.key, names(PLOT.STYLES$fill.colors)
          )
        )
      ) %>%
    arrange(data.type, fill.key) %>%
    droplevels()
  view.lines <- if ("Empirical" %in% active.sources) genome.lines else {
    genome.lines[0, , drop = FALSE]
    }
  view.lines <- view.lines %>%
    filter(as.character(pop) %in% active.populations) %>%
    mutate(pop = factor(as.character(pop), levels = active.populations))
  return(list(points = view.points, genome.lines = view.lines))
  }


# construct one configured selected-chromosome diversity plot
make.diversity.plot <- function(
    points, genome.lines, styles, data.types, tag
  ) {
  view <- filter.diversity.plot.view(
    points, genome.lines, data.types, tag
    )
  points <- view$points
  genome.lines <- view$genome.lines
  points <- filter(points, mask == "Intergenic")
  genome.lines <- filter(genome.lines, mask == "Intergenic")
  dodge <- position_dodge(width = CATEGORICAL.BAR.DODGE)
  fill.keys <- levels(points$fill.key)
  plot <- ggplot(
    points,
    aes(
      x = chrom, y = estimate, fill = fill.key,
      group = interaction(pop, data.type)
      )
    )
    
  plot <- plot +
    geom_col(
      position = dodge, width = CATEGORICAL.BAR.WIDTH,
      color = "black", linewidth = CATEGORICAL.BAR.LINEWIDTH
      ) +
    geom_errorbar(
      data = points,
      aes(ymin = lower, ymax = upper),
      position = dodge, width = 0,
      linewidth = CATEGORICAL.BAR.LINEWIDTH,
      na.rm = TRUE
      ) +
    facet_grid(
      stat ~ ., scales = "free_y",
      labeller = labeller(
        stat = c(pi = "π", theta = "θ[w]")
        )
      ) +
    scale_color_manual(values = styles$empirical.colors) +
    scale_fill_manual(
      values = styles$fill.colors[fill.keys],
      breaks = fill.keys,
      labels = styles$fill.labels[fill.keys],
      limits = fill.keys,
      drop = TRUE
      )
  
  plot <- plot +
    geom_hline(
      data = genome.lines,
      aes(yintercept = estimate, color = pop),
      linetype = "longdash", linewidth = CATEGORICAL.BAR.LINEWIDTH
      ) +
    scale_y_continuous(
      labels = scales::label_number(accuracy = 0.00001)
      )
  
  plot <- plot +
    labs(
      x = "Chromosome", y = NULL,
      title = "Genetic Diversity Across Selected Chromosomes",
      color = NULL, fill = NULL, shape = NULL
      ) +
    guides(
      color = "none",
      fill = guide_legend(order = 1, nrow = 1, byrow = TRUE)
      ) +
    theme(legend.title = element_blank())

  return(apply.standard.plot.theme(plot))
  }


# analysis ----


# read all available chromosome-level simulation diversity estimates
diversity.sim.tc <- read.diversity.chromosomes(
  SIM.TC.DATA.DIR, "pi_theta_stats.chr{chrom}.parquet",
  CHROMOSOMES, "Simulation_2T12Consistent"
  )
diversity.simDown.tc.intergenic <- read.diversity.chromosomes(
  SIMDOWN.TC.DATA.DIR,
  "pi_theta_stats_intergenic.chr{chrom}.parquet",
  CHROMOSOMES, "Simulation_2T12Consistent_simDown", "Intergenic"
  )
diversity.simDown.tc.full.callable <- read.diversity.chromosomes(
  SIMDOWN.TC.DATA.DIR,
  "pi_theta_stats_full_callable_chrom.chr{chrom}.parquet",
  CHROMOSOMES, "Simulation_2T12Consistent_simDown", "Full callable"
  )
diversity.sim.lg <- read.diversity.chromosomes(
  SIM.LG.DATA.DIR, "pi_theta_stats.chr{chrom}.parquet",
  CHROMOSOMES, "Simulation_largeGrowth"
  )
diversity.simDown.lg.intergenic <- read.diversity.chromosomes(
  SIMDOWN.LG.DATA.DIR,
  "pi_theta_stats_intergenic.chr{chrom}.parquet",
  CHROMOSOMES, "Simulation_largeGrowth_simDown", "Intergenic"
  )
diversity.simDown.lg.full.callable <- read.diversity.chromosomes(
  SIMDOWN.LG.DATA.DIR,
  "pi_theta_stats_full_callable_chrom.chr{chrom}.parquet",
  CHROMOSOMES, "Simulation_largeGrowth_simDown", "Full callable"
  )

# read all available chromosome-level empirical diversity estimates
diversity.emp.intergenic.chromosome <- read.diversity.chromosomes(
  EMPIRICAL.DATA.DIR,
  "pi_theta_stats_intergenic.chr{chrom}.parquet",
  CHROMOSOMES, "Empirical", "Intergenic"
  )
diversity.emp.full.callable.chromosome <- read.diversity.chromosomes(
  EMPIRICAL.DATA.DIR,
  "pi_theta_stats_full_callable_chrom.chr{chrom}.parquet",
  CHROMOSOMES, "Empirical", "Full callable"
  )

# read genome-wide empirical diversity references
diversity.emp.intergenic.genome <- read.diversity.genome(
  EMPIRICAL.DATA.DIR, "pi_theta_stats_intergenic.parquet",
  "Intergenic"
  )
diversity.emp.full.callable.genome <- read.diversity.genome(
  EMPIRICAL.DATA.DIR,
  "pi_theta_stats_full_callable_chrom.parquet", "Full callable"
  )

# summarize sources and construct all configured diversity views
diversity.population.contrast.simulation <- bind_rows(
  mutate(diversity.sim.tc, mask = "Intergenic"),
  diversity.simDown.tc.intergenic,
  mutate(diversity.sim.lg, mask = "Intergenic"),
  diversity.simDown.lg.intergenic
  )
diversity.population.contrast.empirical <- bind_rows(
  diversity.emp.intergenic.chromosome, diversity.emp.intergenic.genome
  )
diversity.population.contrast.tables <- make.diversity.population.contrast.tables(
  diversity.population.contrast.simulation, CHROMOSOMES,
  DIVERSITY.POPULATION.CONTRAST.BOOTSTRAP.REPLICATES
  )
diversity.population.contrast.empirical.tables <-
  make.diversity.population.empirical.contrasts(
    diversity.population.contrast.empirical, CHROMOSOMES
    )
diversity.population.contrast.plots <- imap(
  DIVERSITY.POPULATION.CONTRAST.SOURCE.GROUPS,
  function(sources, tag) {
    list(
      `95` = make.diversity.contrast.plot(
        prepare.diversity.population.contrast.plot.data(
          diversity.population.contrast.tables, "95", CHROMOSOMES, sources
          ),
        "95", "population"
        ),
      bonferroni = make.diversity.contrast.plot(
        prepare.diversity.population.contrast.plot.data(
          diversity.population.contrast.tables, "bonferroni", CHROMOSOMES, sources
          ),
        "bonferroni", "population"
        )
      )
    }
  )
diversity.population.contrast.plots$empirical <-
  make.diversity.population.empirical.contrast.plot(
    prepare.diversity.population.empirical.contrast.plot.data(
      diversity.population.contrast.empirical.tables, CHROMOSOMES
      )
    )
diversity.simulation.contrast.tables <- make.diversity.simulation.contrast.tables(
  diversity.population.contrast.simulation, CHROMOSOMES,
  DIVERSITY.POPULATION.CONTRAST.BOOTSTRAP.REPLICATES
  )
diversity.simulation.contrast.plots <- list(
  `95` = make.diversity.contrast.plot(
    prepare.diversity.simulation.contrast.plot.data(
      diversity.simulation.contrast.tables, "95", CHROMOSOMES
      ),
    "95", "simulation"
    ),
  bonferroni = make.diversity.contrast.plot(
    prepare.diversity.simulation.contrast.plot.data(
      diversity.simulation.contrast.tables, "bonferroni", CHROMOSOMES
      ),
    "bonferroni", "simulation"
    )
  )
diversity.simulation.summary <- bind_rows(
  diversity.sim.tc,
  diversity.simDown.tc.intergenic,
  diversity.simDown.tc.full.callable,
  diversity.sim.lg,
  diversity.simDown.lg.intergenic,
  diversity.simDown.lg.full.callable
  ) %>%
  filter(
    !(data.type %in% c(
      "Simulation_largeGrowth", "Simulation_largeGrowth_simDown"
      )) | pop == "ADX"
  ) %>%
  summarize.simulation.diversity()
diversity.plot.data <- build.diversity.plot.data(
  diversity.simulation.summary,
  bind_rows(
    diversity.emp.intergenic.chromosome,
    diversity.emp.full.callable.chromosome
    ),
  bind_rows(
    diversity.emp.intergenic.genome, diversity.emp.full.callable.genome
    ),
  SELECTED.CHROMOSOMES
  )
diversity.bootstrap.plots <- imap(list(
  tc.tcd.1kg = PLOT.CONFIGS$tc.tcd.1kg,
  all.datatypes.adx.asw = PLOT.CONFIGS$all.datatypes.adx.asw
  ), function(data.types, tag) {
  return(make.diversity.plot(
    diversity.plot.data$points,
    diversity.plot.data$genome.lines,
    PLOT.STYLES, data.types,
    tag
    ))
  })

# persist every plot before printing figures at the end of the script
dir.create(OUTPUT.DIR, recursive = TRUE, showWarnings = FALSE)
saveRDS(diversity.bootstrap.plots$tc.tcd.1kg, file.path(
  OUTPUT.DIR, "diversity.bootstrap.tc.tcd.1kg.rds"
  ))
saveRDS(diversity.bootstrap.plots$all.datatypes.adx.asw, file.path(
  OUTPUT.DIR, "diversity.bootstrap.all.datatypes.adx.asw.rds"
  ))
saveRDS(diversity.population.contrast.plots$tc.tcd$`95`, file.path(
  OUTPUT.DIR, "diversity.population.contrasts.tc.tcd.95.rds"
  ))
saveRDS(diversity.population.contrast.plots$tc.tcd$bonferroni, file.path(
  OUTPUT.DIR, "diversity.population.contrasts.tc.tcd.bonferroni.rds"
  ))
saveRDS(diversity.population.contrast.plots$lg.lgd$`95`, file.path(
  OUTPUT.DIR, "diversity.population.contrasts.lg.lgd.95.rds"
  ))
saveRDS(diversity.population.contrast.plots$lg.lgd$bonferroni, file.path(
  OUTPUT.DIR, "diversity.population.contrasts.lg.lgd.bonferroni.rds"
  ))
saveRDS(diversity.population.contrast.plots$empirical, file.path(
  OUTPUT.DIR, "diversity.population.contrasts.empirical.rds"
  ))
saveRDS(diversity.simulation.contrast.plots$`95`, file.path(
  OUTPUT.DIR, "diversity.simulation.contrasts.95.rds"
  ))
saveRDS(diversity.simulation.contrast.plots$bonferroni, file.path(
  OUTPUT.DIR, "diversity.simulation.contrasts.bonferroni.rds"
  ))

print(diversity.bootstrap.plots$tc.tcd.1kg)
print(diversity.bootstrap.plots$all.datatypes.adx.asw)
print(diversity.population.contrast.plots$tc.tcd$`95`)
print(diversity.population.contrast.plots$tc.tcd$bonferroni)
print(diversity.population.contrast.plots$lg.lgd$`95`)
print(diversity.population.contrast.plots$lg.lgd$bonferroni)
print(diversity.population.contrast.plots$empirical)
print(diversity.simulation.contrast.plots$`95`)
print(diversity.simulation.contrast.plots$bonferroni)

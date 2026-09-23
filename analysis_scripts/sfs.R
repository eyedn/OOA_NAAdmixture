# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology
# Mooney Lab
# ___
# sfs.R
# ______________________________________________________________________________

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
SOURCE.DISPLAY.LEVELS <- c("T.C.", "T.C.D.", "L.G.", "L.G.D.", "Emp.")
SOURCE.LABELS <- setNames(SOURCE.DISPLAY.LEVELS, SOURCE.LEVELS)
POPULATION.LEVELS <- c("AFR", "ADX", "EUR", "YRI", "ASW", "CEU")
PLOT.CONFIGS <- list(
  TC.1kG = SOURCE.LEVELS[c(1, 5)],
  TC.TCD = SOURCE.LEVELS[c(1, 2)],
  TCD.1kG = SOURCE.LEVELS[c(2, 5)],
  onlyADX = SOURCE.LEVELS[1:4],
  tcd.1kg = SOURCE.LEVELS[c(2, 5)],
  tc.tcd.1kg = SOURCE.LEVELS[c(1, 2, 5)],
  all.datatypes.adx.asw = SOURCE.LEVELS
  )
RANDOM.SEED <- 123L
BOOTSTRAP.REPLICATES <- 1000L
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
CATEGORICAL.BAR.DODGE <- 0.9
CATEGORICAL.BAR.WIDTH <- 0.8
CATEGORICAL.BAR.LINEWIDTH <- 1
DENSE.BAR.WIDTH.MULTIPLIER <- 0.8
DENSE.BAR.LINEWIDTH <- 0.75
PLOT.STYLES <- list(series.labels = SOURCE.LABELS)


# internal functions ----


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


# return source-specific population pairs used in SFS diagnostics
population.pair.config <- function() {
  config <- tribble(
    ~data.type, ~comparison, ~pair.left, ~pair.right,
    "Simulation_2T12Consistent", "AFR - EUR", "AFR", "EUR",
    "Simulation_2T12Consistent_simDown", "AFR - EUR", "AFR", "EUR",
    "Empirical", "YRI - CEU", "YRI", "CEU"
    )
  return(config)
  }


# prepare singleton denominators for simulation and empirical population pairs
prepare.singleton.diagnostics <- function(data) {
  pair.config <- population.pair.config()
  diagnostic.data <- data %>%
    mutate(
      data.type = as.character(data.type),
      chrom = as.character(chrom)
      ) %>%
    inner_join(pair.config, by = "data.type") %>%
    filter(pop %in% c(pair.left, pair.right)) %>%
    mutate(
      bin.range = case_when(
        minor.allele.count == 1L ~ "Singletons (bin 1)",
        minor.allele.count <= DISPLAY.BIN.MAX ~ "Bins 2-15",
        TRUE ~ "Bins 16-50"
        ),
      bin.range = factor(
        bin.range,
        levels = c("Singletons (bin 1)", "Bins 2-15", "Bins 16-50")
        )
      )
  group.columns <- c(
    "data.set", "data.type", "rep", "chrom", "comparison", "pair.left",
    "pair.right", "pop", "series"
    )
  composition <- diagnostic.data %>%
    group_by(across(all_of(c(group.columns, "bin.range")))) %>%
    summarize(count = sum(count), .groups = "drop")
  partition.check <- composition %>%
    group_by(across(all_of(group.columns))) %>%
    summarize(
      n.partitions = n(),
      partition.total = sum(count),
      .groups = "drop"
      )
  segregating.totals <- diagnostic.data %>%
    group_by(across(all_of(group.columns))) %>%
    summarize(segregating.total = sum(count), .groups = "drop")
  if (any(partition.check$n.partitions != 3L)) {
    stop("Every singleton diagnostic spectrum must have three partitions")
    }
  if (!isTRUE(all.equal(
    partition.check$partition.total,
    segregating.totals$segregating.total
    ))) {
    stop("Singleton diagnostic partitions do not sum to bins 1-50")
    }
  composition <- composition %>%
    left_join(partition.check, by = group.columns) %>%
    left_join(segregating.totals, by = group.columns)
  primary.singletons <- diagnostic.data %>%
    filter(minor.allele.count == 1L) %>%
    select(
      all_of(group.columns),
      primary.singleton.proportion = proportion
      )
  paired <- composition %>%
    filter(bin.range == "Singletons (bin 1)") %>%
    transmute(
      across(all_of(group.columns)),
      singleton.proportion = count / segregating.total
      ) %>%
    left_join(primary.singletons, by = group.columns)
  if (!isTRUE(all.equal(
    paired$singleton.proportion,
    paired$primary.singleton.proportion
    ))) {
    stop("Singleton diagnostic proportions do not match the primary SFS")
    }
  paired.populations <- paired %>%
    group_by(data.type, rep, chrom, comparison, pair.left, pair.right) %>%
    summarize(
      valid = identical(
        sort(pop),
        sort(c(first(pair.left), first(pair.right)))
        ),
      .groups = "drop"
      )
  if (!all(paired.populations$valid)) {
    stop("Every singleton diagnostic needs its configured population pair")
    }
  simulation.sources <- c(
    "Simulation_2T12Consistent",
    "Simulation_2T12Consistent_simDown"
    )
  diagnostics <- list(
    simulation = list(
      composition = composition %>% filter(data.type %in% simulation.sources),
      paired = paired %>% filter(data.type %in% simulation.sources)
      ),
    empirical = list(
      composition = composition %>% filter(data.type == "Empirical"),
      paired = paired %>% filter(data.type == "Empirical")
      )
    )
  return(diagnostics)
  }


# calculate per-replicate population differences over displayed SFS bins
prepare.population.differences <- function(data) {
  pair.config <- population.pair.config()
  pair.columns <- c(
    "data.set", "data.type", "rep", "chrom", "comparison", "pair.left",
    "pair.right", "minor.allele.count"
    )
  paired.data <- data %>%
    mutate(
      data.type = as.character(data.type),
      chrom = as.character(chrom)
      ) %>%
    inner_join(pair.config, by = "data.type") %>%
    filter(
      pop %in% c(pair.left, pair.right),
      minor.allele.count <= DISPLAY.BIN.MAX
      ) %>%
    mutate(pair.member = if_else(pop == pair.left, "left", "right"))
  pair.check <- paired.data %>%
    group_by(across(all_of(pair.columns))) %>%
    summarize(n.populations = n_distinct(pair.member), .groups = "drop")
  if (any(pair.check$n.populations != 2L)) {
    stop("Every population difference needs both configured populations")
    }
  differences <- paired.data %>%
    select(all_of(pair.columns), pair.member, count, proportion) %>%
    pivot_wider(
      names_from = pair.member,
      values_from = c(count, proportion)
      ) %>%
    transmute(
      across(all_of(pair.columns)),
      count.difference = count_left - count_right,
      proportion.difference = proportion_left - proportion_right
      ) %>%
    pivot_longer(
      cols = c(count.difference, proportion.difference),
      names_to = "measure",
      values_to = "difference"
      ) %>%
    mutate(
      comparison = factor(
        comparison,
        levels = c("AFR - EUR", "YRI - CEU")
        ),
      measure = recode(
        measure,
        count.difference = "Count difference",
        proportion.difference = "Proportion difference"
        ),
      measure = factor(
        measure,
        levels = c("Count difference", "Proportion difference")
        ),
      data.type = factor(data.type, levels = SOURCE.LEVELS)
      )
  return(differences)
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
    filter(tag != "onlyADX" | pop == "ADX") %>%
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
      width = 0.25, linewidth = DENSE.BAR.LINEWIDTH,
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
    guides(fill = guide_legend(order = 1, nrow = 1, byrow = TRUE)) +
    theme_bw(base_size = PLOT.BASE.SIZE) +
    theme(
      legend.position = "top", legend.direction = "horizontal", 
      legend.box = "horizontal", panel.grid.minor = element_blank()
    )
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
  return(plot)
  }


# build simulation count compositions with empirical population references
make.singleton.composition.plot <- function(data, empirical) {
  source.labels <- PLOT.STYLES$series.labels
  data <- data %>%
    mutate(
      data.type = factor(
        as.character(data.type),
        levels = order.active.levels(data.type, SOURCE.LEVELS)
        ),
      pop = factor(
        as.character(pop),
        levels = order.active.levels(pop, POPULATION.LEVELS)
        )
      )
  empirical <- empirical %>%
    mutate(pop = factor(
      as.character(pop),
      levels = order.active.levels(pop, POPULATION.LEVELS)
      ))
  plot <- ggplot(
    data,
    aes(
      x = factor(rep),
      y = count,
      fill = pop,
      group = pop
      )
    ) +
    geom_col(position = SFS.DODGE, width = 0.8) +
    geom_hline(
      data = empirical %>% select(bin.range, pop, count),
      aes(yintercept = count, color = pop),
      linewidth = 0.8, linetype = "dashed"
      ) +
    facet_grid(
      rows = vars(data.type),
      cols = vars(bin.range),
      labeller = labeller(data.type = source.labels),
      drop = FALSE
      ) +
    scale_fill_manual(values = c(
      AFR = "#56B4E9",
      EUR = "#FB8072"
      )) +
    scale_color_manual(values = c(
      YRI = "#EEC4DC",
      CEU = "#BB437E"
      )) +
    labs(
      x = "Replicate",
      y = "Segregating-site count",
      fill = "Simulation population",
      color = "Empirical reference",
      title = "Singleton denominator composition"
      ) +
    theme_bw(base_size = 18) +
    theme(
      legend.position = "top",
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1)
      )
  return(plot)
  }


# build a two-row SFS difference plot with all source pairs overlaid
make.population.difference.plot <- function(data) {
  source.labels <- PLOT.STYLES$series.labels
  source.colors <- c(
    Simulation_2T12Consistent = "#56B4E9",
    Simulation_2T12Consistent_simDown = "#6F55B5",
    Empirical = "#E44B8D"
    )
  data <- data %>%
    mutate(data.type = factor(
      as.character(data.type),
      levels = order.active.levels(data.type, SOURCE.LEVELS)
      ))
  plot <- ggplot(
    data,
    aes(
      x = minor.allele.count,
      y = difference,
      color = data.type,
      group = interaction(data.type, rep, chrom)
      )
    ) +
    geom_hline(yintercept = 0, color = "grey60", linewidth = 0.5) +
    geom_line(alpha = 0.35) +
    stat_summary(
      aes(group = data.type),
      fun = mean,
      geom = "line",
      linewidth = 1
      ) +
    facet_grid(measure ~ ., scales = "free_y", drop = FALSE) +
    scale_x_continuous(breaks = seq_len(DISPLAY.BIN.MAX)) +
    scale_color_manual(
      values = source.colors,
      breaks = levels(data$data.type),
      labels = source.labels[levels(data$data.type)]
      ) +
    labs(
      x = "Minor allele count",
      y = "Left population minus right population",
      color = NULL,
      title = "Population differences across shown SFS bins",
      subtitle = paste(
        "Thin lines are replicates; thick lines are source means.",
        "Simulations: AFR - EUR. Empirical: YRI - CEU."
        )
      ) +
    theme_bw(base_size = 18) +
    theme(
      legend.position = "top",
      panel.grid.minor = element_blank()
      )
  return(plot)
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
# singleton.diagnostics <- prepare.singleton.diagnostics(sfs.data)
# population.differences <- prepare.population.differences(sfs.data)

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
# singleton.composition.plot <- make.singleton.composition.plot(
#   singleton.diagnostics$simulation$composition,
#   singleton.diagnostics$empirical$composition
#   )
# population.difference.plot <- make.population.difference.plot(
#   population.differences
#   )

# Legacy diagnostics remain available above but are not emitted in this refresh.
# Persist chromosome-1 bootstrap count and proportion plots before printing.
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

print(sfs.bootstrap.count.plots$tcd.1kg)
print(sfs.bootstrap.count.plots$all.datatypes.adx.asw)
print(sfs.bootstrap.proportion.plots$tcd.1kg)
print(sfs.bootstrap.proportion.plots$all.datatypes.adx.asw)

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
library(furrr)
library(nanoparquet)
library(scales)


SIM.SMALL.DATA.DIR <- "~/scratch/OOA_NAAdmixture_small/stats"
SIMDOWN.SMALL.DATA.DIR <- "~/scratch/OOA_NAAdmixture_smallOnekgDownsample/stats"
SIM.LARGE.DATA.DIR <- "~/scratch/OOA_NAAdmixture_large/stats"
SIMDOWN.LARGE.DATA.DIR <- "~/scratch/OOA_NAAdmixture_largeOnekgDownsample/stats"
EMPIRICAL.DATA.DIR <- "~/scratch/OOA_NAAdmixture_1kG/stats"
# CHROMOSOMES <- as.character(1:22)
CHROMOSOMES <- c("1")
DISPLAY.BIN.MAX <- 15
SFS.PROJECTION.ALLELE.COUNT <- 100
SOURCE.LEVELS <- c(
  "Simulation_small", "Simulation_large", "Simulation_small_simDown",
  "Simulation_large_simDown", "Empirical"
)
SFS.FACET.LEVELS <- c("1", "all")
SFS.SERIES.LEVELS <- c(
  "small AFR", "small ADX", "small EUR",
  "large ADX", "small simDown AFR", "small simDown ADX",
  "small simDown EUR", "large simDown ADX",
  "empirical YRI", "empirical ASW", "empirical CEU"
)
SFS.COLORS <- c(
  "small AFR" = "#9BD5F2",
  "small ADX" = "#9A83CE",
  "small EUR" = "#FBB4AE",
  "small simDown AFR" = "#56B4E9",
  "small simDown ADX" = "#6F55B5",
  "small simDown EUR" = "#FB8072",
  "large ADX" = "#32146F",
  "large simDown ADX" = "#4B1FA8",
  "empirical YRI" = "#EEC4DC",
  "empirical ASW" = "#E44B8D",
  "empirical CEU" = "#BB437E"
)
SFS.SERIES.LABELS <- c(
  "small AFR" = "Sm. Sim. AFR",
  "small ADX" = "Sm. Sim. ADX",
  "small EUR" = "Sm. Sim. EUR",
  "small simDown AFR" = "Sm. D. Sim. AFR",
  "small simDown ADX" = "Sm. D. Sim. ADX",
  "small simDown EUR" = "Sm. D. Sim. EUR",
  "large ADX" = "Lg. Sim. ADX",
  "large simDown ADX" = "Lg. D. Sim. ADX",
  "empirical YRI" = "Emp. YRI",
  "empirical ASW" = "Emp. ASW",
  "empirical CEU" = "Emp. CEU"
)
SFS.DODGE <- position_dodge(width = 0.9)
PLOT.STYLES <- list(series.labels = c(
  Simulation_small = "Sm. Sim.",
  Simulation_large = "Lg. Sim.",
  Simulation_small_simDown = "Sm. D. Sim.",
  Simulation_large_simDown = "Lg. D. Sim.",
  Empirical = "Emp."
))


# internal functions ----


# describe the shared projection and simulation uncertainty concisely
sfs.plot.subtitle <- function() {
  subtitle <- paste0(
    "Projected to ", SFS.PROJECTION.ALLELE.COUNT,
    " alleles · error bars show ±2 SD"
  )
  return(subtitle)
}

# retain source-specific populations before SFS normalization
apply.sfs.source.contract <- function(data) {
  retained <- data %>%
    filter(
      (data.type %in% c(
        "Simulation_small", "Simulation_small_simDown"
      ) & pop %in% c("AFR", "ADX", "EUR")) |
        (data.type %in% c(
          "Simulation_large", "Simulation_large_simDown"
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
        "Simulation_large", "Simulation_large_simDown"
      ),
      "large", "small"
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
    mutate(data.type = if_else(data.set == "small", "Simulation_small",
                               "Simulation_large")) %>%
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
        data.set = if_else(grepl("large", data.type), "large", "small"),
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
    filter(chrom %in% SFS.FACET.LEVELS) %>%
    mutate(
      series = case_when(
        data.type == "Simulation_small" ~ paste("small", pop),
        data.type == "Simulation_small_simDown" ~
          paste("small simDown", pop),
        data.type == "Simulation_large" ~ paste("large", pop),
        data.type == "Simulation_large_simDown" ~
          paste("large simDown", pop),
        data.type == "Empirical" ~ paste("empirical", pop)
      ),
      chrom = factor(chrom, levels = SFS.FACET.LEVELS),
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
  summary <- data %>%
    group_by(
      data.set, data.type, chrom, pop, series, minor.allele.count
    ) %>%
    summarize(
      mean = mean(.data[[value.column]]),
      sd = if_else(
        first(data.set) == "empirical",
        NA_real_,
        sd(.data[[value.column]])
      ),
      .groups = "drop"
    ) %>%
    mutate(
      lower = pmax(0, mean - 2 * sd),
      upper = mean + 2 * sd,
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


# build one shared dodged-bar SFS plot
filter.plot.view <- function(data, view) {
  sources <- switch(
    view,
    small_empirical = c(
      "Simulation_small", "Simulation_small_simDown", "Empirical"
    ),
    simulated = c(
      "Simulation_small", "Simulation_large", "Simulation_small_simDown",
      "Simulation_large_simDown"
    ),
    stop("Unsupported SFS plot view: ", view)
  )
  filtered <- data %>%
    filter(as.character(data.type) %in% sources) %>%
    filter(view != "simulated" | pop == "ADX") %>%
    mutate(
      data.type = factor(as.character(data.type), levels = sources),
      series = factor(
        as.character(series),
        levels = SFS.SERIES.LEVELS[SFS.SERIES.LEVELS %in% series]
      ),
      chrom = factor(
        as.character(chrom),
        levels = SFS.FACET.LEVELS[SFS.FACET.LEVELS %in% chrom]
      )
    ) %>%
    arrange(data.type, series) %>%
    droplevels()
  return(filtered)
}


# build one scoped shared dodged-bar SFS plot
make.sfs.plot <- function(data, value.column, y.label, pseudo.log, view) {
  displayed <- filter.plot.view(data, view) %>%
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
    geom_col(position = SFS.DODGE, width = 0.9) +
    geom_errorbar(
      aes(ymin = lower, ymax = upper),
      position = SFS.DODGE,
      width = 0.25, linewidth = 0.5,
      na.rm = TRUE
    ) +
    facet_wrap(~chrom, nrow = 1, drop = FALSE) +
    scale_x_continuous(
      breaks = seq_len(DISPLAY.BIN.MAX),
      limits = c(0.5, DISPLAY.BIN.MAX + 0.5)
    ) +
    scale_fill_manual(
      values = SFS.COLORS[series.keys],
      breaks = series.keys,
      labels = SFS.SERIES.LABELS[series.keys],
      limits = series.keys,
      drop = FALSE
    ) +
    labs(
      x = "Minor allele count",
      y = y.label,
      subtitle = sfs.plot.subtitle(),
      fill = NULL
    ) +
    theme_bw(base_size = 24) +
    theme(
      legend.position = "top", legend.direction = "horizontal", 
      legend.box = "horizontal", panel.grid.minor = element_blank()
    )
  if (pseudo.log) {
    plot <- plot + scale_y_continuous(trans = pseudo_log_trans())
  } else {
    plot <- plot + scale_y_continuous()
  }
  return(plot)
}


# read standardized projected chromosome producer outputs
read.sfs.inputs <- function(
    sim.small.dir, simDown.small.dir, sim.large.dir,
    simDown.large.dir, empirical.dir
) {
  read.chromosomes <- function(directory) {
    previous.plan <- future::plan()
    on.exit(future::plan(previous.plan), add = TRUE)
    future::plan(future::multisession)
    tables <- future_map_dfr(CHROMOSOMES, function(chrom) {
      table <- read_parquet(
        file.path(
          path.expand(directory),
          paste0("sfs.chr", chrom, ".parquet")
        )
      )
      table$chrom <- chrom
      return(table)
    }, .options = furrr_options(seed = TRUE))
    return(tables)
  }
  simulation <- bind_rows(
    read.chromosomes(sim.small.dir) %>% mutate(data.set = "small"),
    read.chromosomes(sim.large.dir) %>% mutate(data.set = "large")
  )
  simDown <- bind_rows(
    read.chromosomes(simDown.small.dir) %>%
      mutate(data.type = "Simulation_small_simDown"),
    read.chromosomes(simDown.large.dir) %>%
      mutate(data.type = "Simulation_large_simDown")
  )
  empirical <- read.chromosomes(empirical.dir)
  return(list(
    simulation = simulation,
    simDown = simDown,
    empirical = empirical
  ))
}


# analysis ----


sfs.inputs <- read.sfs.inputs(
  SIM.SMALL.DATA.DIR,
  SIMDOWN.SMALL.DATA.DIR,
  SIM.LARGE.DATA.DIR,
  SIMDOWN.LARGE.DATA.DIR,
  EMPIRICAL.DATA.DIR
)
sfs.data <- prepare.sfs.analysis(
  sfs.inputs$simulation,
  sfs.inputs$simDown,
  sfs.inputs$empirical
)
sfs.summaries <- summarize.sfs.analysis(sfs.data)

sfs.small.empirical.count.plot <- make.sfs.plot(
  sfs.summaries$count,
  "mean",
  "Projected site count",
  TRUE, "small_empirical"
)
sfs.simulated.count.plot <- make.sfs.plot(
  sfs.summaries$count,
  "mean",
  "Projected site count",
  TRUE, "simulated"
)
sfs.small.empirical.proportion.plot <- make.sfs.plot(
  sfs.summaries$proportion,
  "mean",
  "Proportion of segregating sites",
  FALSE, "small_empirical"
)
sfs.simulated.proportion.plot <- make.sfs.plot(
  sfs.summaries$proportion,
  "mean",
  "Proportion of segregating sites",
  FALSE, "simulated"
)

# print(sfs.small.empirical.count.plot)
# print(sfs.simulated.count.plot)
print(sfs.small.empirical.proportion.plot)
print(sfs.simulated.proportion.plot)

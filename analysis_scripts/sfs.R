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
# Reason: self-contained helpers accompany sequential I/O orchestration.

# set up ----
library(tidyverse)
library(furrr)
library(nanoparquet)
library(scales)


SIM.SMALL.DATA.DIR <- "~/scratch/OOA_NAAdmixture_small/stats"
SIM.LARGE.DATA.DIR <- "~/scratch/OOA_NAAdmixture_large/stats"
EMPIRICAL.DATA.DIR <- "~/scratch/OOA_NAAdmixture_1kG/stats"
CHROMOSOMES <- as.character(1:22)
DISPLAY.BIN.MAX <- 10
SFS.PROJECTION.ALLELE.COUNT <- 100
SFS.FACET.LEVELS <- c("all", "1", "5", "10", "14", "18", "22")
SFS.SERIES.LEVELS <- c(
  "small AFR", "small ADX", "small EUR",
  "large AFR", "large ADX", "large EUR",
  "empirical YRI", "empirical ASW", "empirical CEU"
)
SFS.COLORS <- c(
  "small AFR" = "#9BD5F2",
  "small ADX" = "#9A83CE",
  "small EUR" = "#FBB4AE",
  "large AFR" = "#0072B2",
  "large ADX" = "#32146F",
  "large EUR" = "#D94A3A",
  "empirical YRI" = "#EEC4DC",
  "empirical ASW" = "#E44B8D",
  "empirical CEU" = "#BB437E"
)
SFS.DODGE <- position_dodge(width = 0.9)


# internal functions ----


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
    group_by(data.set, rep, chrom, pop, minor.allele.count) %>%
    summarize(count = sum(count), .groups = "drop")
  return(folded)
}


# add a genome spectrum by summing projected chromosome spectra
add.simulation.genome <- function(data) {
  genome <- data %>%
    group_by(data.set, rep, pop, derived_allele_count) %>%
    summarize(count = sum(count), .groups = "drop") %>%
    mutate(chrom = "all", .before = pop)
  return(bind_rows(data, genome))
}


# add a genome spectrum by summing projected empirical chromosomes
add.empirical.genome <- function(data) {
  genome <- data %>%
    group_by(data.set, rep, pop, minor.allele.count) %>%
    summarize(count = sum(count), .groups = "drop") %>%
    mutate(chrom = "all", .before = pop)
  return(bind_rows(data, genome))
}


# standardize producers, fold simulation once, and normalize full spectra
prepare.sfs.analysis <- function(simulation, empirical) {
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
  if (any(empirical$projection_allele_count !=
          SFS.PROJECTION.ALLELE.COUNT)) {
    stop("Empirical SFS projection metadata do not match the analysis")
  }

  simulation <- simulation %>%
    mutate(chrom = as.character(chrom)) %>%
    validate.complete.sfs(
      "derived_allele_count",
      0:SFS.PROJECTION.ALLELE.COUNT
    ) %>%
    add.simulation.genome() %>%
    fold.simulation.sfs()
  empirical <- empirical %>%
    mutate(data.set = "empirical", chrom = as.character(chrom)) %>%
    validate.complete.sfs(
      "minor_allele_count",
      0:(SFS.PROJECTION.ALLELE.COUNT / 2)
    ) %>%
    rename(minor.allele.count = minor_allele_count) %>%
    add.empirical.genome() %>%
    select(data.set, rep, chrom, pop, minor.allele.count, count)

  prepared <- bind_rows(simulation, empirical) %>%
    filter(chrom %in% SFS.FACET.LEVELS) %>%
    mutate(
      series = paste(data.set, pop),
      chrom = factor(chrom, levels = SFS.FACET.LEVELS),
      series = factor(series, levels = SFS.SERIES.LEVELS)
    )
  if (any(is.na(prepared$series))) {
    stop("SFS inputs contain an unsupported data-set/population series")
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
  return(prepared)
}


# summarize simulation replicates and retain empirical NA intervals
summarize.one.sfs.value <- function(data, value.column) {
  summary <- data %>%
    group_by(data.set, chrom, series, minor.allele.count) %>%
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
make.sfs.plot <- function(data, value.column, y.label, pseudo.log) {
  displayed <- data %>%
    filter(minor.allele.count <= DISPLAY.BIN.MAX)
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
      width = 0.25,
      na.rm = TRUE
    ) +
    facet_wrap(~chrom, nrow = 1, drop = FALSE) +
    scale_x_continuous(
      breaks = seq_len(DISPLAY.BIN.MAX),
      limits = c(0.5, DISPLAY.BIN.MAX + 0.5)
    ) +
    scale_fill_manual(values = SFS.COLORS, drop = FALSE) +
    labs(
      x = "Minor allele count",
      y = y.label,
      fill = NULL
    ) +
    theme_bw(base_size = 18) +
    theme(
      legend.position = "top",
      axis.text.x = element_text(angle = 90, vjust = 0.5)
    )
  if (pseudo.log) {
    plot <- plot + scale_y_continuous(trans = pseudo_log_trans())
  } else {
    plot <- plot + scale_y_continuous()
  }
  return(plot)
}


# read standardized projected chromosome producer outputs
read.sfs.inputs <- function(sim.small.dir, sim.large.dir, empirical.dir) {
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
  empirical <- read.chromosomes(empirical.dir)
  return(list(simulation = simulation, empirical = empirical))
}


# analysis ----


sfs.inputs <- read.sfs.inputs(
  SIM.SMALL.DATA.DIR,
  SIM.LARGE.DATA.DIR,
  EMPIRICAL.DATA.DIR
)
sfs.data <- prepare.sfs.analysis(
  sfs.inputs$simulation,
  sfs.inputs$empirical
)
sfs.summaries <- summarize.sfs.analysis(sfs.data)

sfs.count.plot <- make.sfs.plot(
  sfs.summaries$count,
  "mean",
  "Projected site count",
  TRUE
)
sfs.proportion.plot <- make.sfs.plot(
  sfs.summaries$proportion,
  "mean",
  "Proportion of segregating sites",
  FALSE
)

print(sfs.count.plot)
print(sfs.proportion.plot)

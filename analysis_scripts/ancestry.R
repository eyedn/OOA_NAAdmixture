# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology
# Mooney Lab
# ___
# ancestry.R
# ______________________________________________________________________________


# set up ----
setwd("~/OOA_NAAdmixture/analysis_scripts/")
options(scipen = 999)
library(tidyverse)
library(ggh4x)
library(nanoparquet)
source(if (file.exists("analysis_scripts/bootstrap_parallel.R")) {
  "analysis_scripts/bootstrap_parallel.R"
  } else {
  "bootstrap_parallel.R"
  })


ANCESTRY.SIM.TC.DATA.DIR <- "~/scratch/OOA_NAAdmixture_2T12Consistent/stats"
ANCESTRY.SIMDOWN.TC.DATA.DIR <-
  "~/scratch/OOA_NAAdmixture_2T12ConsistentOnekgDownsample/stats"
ANCESTRY.SIM.LG.DATA.DIR <- "~/scratch/OOA_NAAdmixture_largeGrowth/stats"
ANCESTRY.SIMDOWN.LG.DATA.DIR <-
  "~/scratch/OOA_NAAdmixture_largeGrowthOnekgDownsample/stats"
ANCESTRY.EMPIRICAL.DATA.DIR <- "~/scratch/OOA_NAAdmixture_1kG/stats"
ANCESTRY.OUTPUT.DIR <- "/home1/karatas/proj/OOA_NAAdmixture_data"
ANCESTRY.CHROMOSOMES <- as.character(1:22)
ANCESTRY.SELECTED.CHROMOSOMES <- c("1", "10", "20")
ANCESTRY.SIMULATION.K <- 2L
ANCESTRY.EMPIRICAL.K <- 2L
ANCESTRY.RANDOM.SEED <- 123L
ANCESTRY.DOWNSAMPLE.SIZE <- 50L
ANCESTRY.BOOTSTRAP.REPLICATES <- 100000L
# per-statistic Bonferroni family for chromosome-level ancestry comparisons.
ANCESTRY.BOOTSTRAP.CHROMOSOME.FAMILY.SIZE <- 154L
# per-statistic Bonferroni family for chromosome versus genome-wide ASW.
ANCESTRY.BOOTSTRAP.GENOME.ASW.FAMILY.SIZE <- 88L
ANCESTRY.HISTOGRAM.BREAKS <- seq(0, 1, by = 0.05)
ANCESTRY.PLOT.EMPIRICAL.METHOD <- "ADMIXTURE"
ANCESTRY.PLOT.SAMPLE.SET <- "downsampled"
ANCESTRY.PLOT.BASE.SIZE <- 24
ANCESTRY.CATEGORICAL.BAR.DODGE <- 0.9
ANCESTRY.CATEGORICAL.BAR.WIDTH <- 0.8
ANCESTRY.CATEGORICAL.BAR.LINEWIDTH <- 1
ANCESTRY.DENSE.BAR.WIDTH.MULTIPLIER <- 0.8
ANCESTRY.DENSE.BAR.LINEWIDTH <- 0.75
ANCESTRY.SOURCE.LEVELS <- c(
  "Simulation_2T12Consistent", "Simulation_2T12Consistent_simDown",
  "Simulation_largeGrowth", "Simulation_largeGrowth_simDown", "Empirical"
  )
ANCESTRY.SOURCE.LABELS <- setNames(
  c("T.C. ADX", "T.C.D. ADX", "L.G. ADX", "L.G.D. ADX", "ASW"), ANCESTRY.SOURCE.LEVELS
  )
ANCESTRY.PLOT.STYLES <- list(
  source.levels = ANCESTRY.SOURCE.LEVELS,
  base.size = ANCESTRY.PLOT.BASE.SIZE,
  categorical.bar.dodge = ANCESTRY.CATEGORICAL.BAR.DODGE,
  categorical.bar.width = ANCESTRY.CATEGORICAL.BAR.WIDTH,
  categorical.bar.linewidth = ANCESTRY.CATEGORICAL.BAR.LINEWIDTH,
  dense.bar.width.multiplier = ANCESTRY.DENSE.BAR.WIDTH.MULTIPLIER,
  dense.bar.linewidth = ANCESTRY.DENSE.BAR.LINEWIDTH,
  empirical.method = ANCESTRY.PLOT.EMPIRICAL.METHOD,
  sample.set = ANCESTRY.PLOT.SAMPLE.SET,
  colors = c(
    Simulation_2T12Consistent = "#9A83CE",
    Simulation_2T12Consistent_simDown = "#6F55B5",
    Simulation_largeGrowth = "#32146F",
    Simulation_largeGrowth_simDown = "#4B1FA8", Empirical = "#E44B8D"
    ),
  labels = ANCESTRY.SOURCE.LABELS,
  empirical.colors = c(ADMIXTURE = "#E44B8D", fastStructure = "#E44B8D"),
  contrast.colors = c(
    `TC-Emp` = "#9A83CE", `TCD-Emp` = "#6F55B5",
    `LG-Emp` = "#32146F", `LGD-Emp` = "#4B1FA8",
    `TC-TCD` = "#7FC2C5", `LG-LGD` = "#2C777C",
    `TC-LG` = "#002526" 
    ),
  contrast.labels = c(
    `TC-Emp` = "T.C. ADX - ASW", `TCD-Emp` = "T.C.D. ADX- ASW",
    `LG-Emp` = "L.G. ADX - ASW", `LGD-Emp` = "L.G.D. ADX - ASW",
    `TC-TCD` = "T.C. ADX - T.C.D. ADX", `LG-LGD` = "L.G. ADX - L.G.D. ADX",
    `TC-LG` = "T.C. ADX - L.G. ADX"
    )
  )
ANCESTRY.BOOTSTRAP.CONTRASTS <- tribble(
  ~contrast, ~left.source, ~right.source, ~paired,
  "TC-TCD", "TC", "TCD", TRUE,
  "LG-LGD", "LG", "LGD", TRUE,
  "TC-LG", "TC", "LG", FALSE,
  "TC-Emp", "TC", "Emp", FALSE,
  "LG-Emp", "LG", "Emp", FALSE,
  "TCD-Emp", "TCD", "Emp", FALSE,
  "LGD-Emp", "LGD", "Emp", FALSE
  )


# internal functions ----


# apply the common visual treatment to a completed plot.
apply.standard.plot.theme <- function(plot, legend.position = "top") {
  return(plot + theme_bw(base_size = ANCESTRY.PLOT.BASE.SIZE) + theme(
    legend.position = legend.position,
    legend.direction = "horizontal",
    legend.box = "horizontal",
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "plain")
    ))
  }


# return the configured empirical-inference subtitle.
ancestry.inference.subtitle <- function() {
  return(paste("Empirical inference:", ANCESTRY.PLOT.EMPIRICAL.METHOD))
  }


# return a method-tagged output filename.
ancestry.output.filename <- function(stem, extension) {
  return(paste0(stem, ".", tolower(ANCESTRY.PLOT.EMPIRICAL.METHOD), ".", extension))
  }


# return active sources in their canonical display order.
order.active.levels <- function(values, canonical.levels) {
  return(canonical.levels[canonical.levels %in% as.character(values)])
  }


# return the configured chromosome-level inference filename family.
ancestry.inference.file.family <- function(method) {
  if (!method %in% c("ADMIXTURE", "fastStructure")) {
    stop("Unsupported ancestry inference method: ", method)
    }
  return(paste0("ancestry_", method, "_multik.chr{chrom}.parquet"))
  }


# retain the sole adx simulation and asw empirical populations.
apply.ancestry.source.contract <- function(data) {
  retained <- data %>%
    filter(
      (data.type %in% ANCESTRY.SOURCE.LEVELS[1:4] & pop == "ADX") |
        (data.type == "Empirical" & pop == "ASW")
      ) %>%
    mutate(data.type = factor(data.type, levels = ANCESTRY.SOURCE.LEVELS))
  return(retained)
  }


# standardize a source table to the shared ancestry schema.
normalize.ancestry.table <- function(
    data, data.type.input, method, simulation.source.input, k,
    sample.id.column = "sample_id"
  ) {
  data$chrom <- as.character(data$chrom)
  data$sample_id <- as.character(data[[sample.id.column]])
  data$rep <- if ("rep" %in% names(data)) as.numeric(data$rep) else 0
  data$data.type <- data.type.input
  data$method <- method
  data$simulation.source <- simulation.source.input
  if (!"pop" %in% names(data)) data$pop <- NA_character_
  if (!"role" %in% names(data)) data$role <- data$pop
  if (method == "tspop") {
    data$component_1_q <- as.numeric(data$afr_tspop)
    data$component_2_q <- as.numeric(data$eur_tspop)
    data <- select(data, -afr_tspop, -eur_tspop)
    data$k <- 0L
    } else {
    if (!"k" %in% names(data)) data$k <- k
    data <- filter(data, .data$k == .env$k)
    }
  data <- data %>% mutate(
    component_1_q = as.numeric(component_1_q),
    component_2_q = as.numeric(component_2_q)
    )
  return(data)
  }


# read available chromosome files and an optional whole-genome empirical file.
read.ancestry.family <- function(
    data.directory, file.family, chromosomes, data.type.input, method, k,
    simulation.source.input, include.genome = FALSE, genome.file.family = NULL
  ) {
  paths <- file.path(
    path.expand(data.directory),
    vapply(
      chromosomes,
      function(chromosome) gsub("\\{chrom\\}", chromosome, file.family),
      character(1)
      )
    )
  missing <- chromosomes != "all" & !file.exists(paths)
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
  chromosomes <- chromosomes[!missing]
  paths <- paths[!missing]
  if (include.genome) {
    chromosomes <- c(chromosomes, "all")
    paths <- c(
      paths,
      file.path(path.expand(data.directory), genome.file.family)
      )
    }
  data <- map2_dfr(paths, chromosomes, function(path, chrom) {
    table <- read_parquet(path)
    table$chrom <- chrom
    return(normalize.ancestry.table(
      table, data.type.input, method, simulation.source.input, k
      ))
    })
  return(data)
  }


# orient k = 2 components so afr.q consistently represents African ancestry.
orient.ancestry.components <- function(data, grouping.columns) {
  component.map <- data %>%
    group_by(across(all_of(grouping.columns))) %>%
    summarise(
      component.1 = mean(component_1_q, na.rm = TRUE),
      component.2 = mean(component_2_q, na.rm = TRUE), .groups = "drop"
      ) %>%
    mutate(afr.component = if_else(
      component.1 >= component.2, "component_1_q", "component_2_q"
      ))
  oriented <- data %>%
    left_join(component.map, by = grouping.columns) %>%
    mutate(afr.q = if_else(
      afr.component == "component_1_q", component_1_q, component_2_q
      )) %>%
    select(-component.1, -component.2, -afr.component)
  return(oriented)
  }


# select reproducible simulation ids within each source, replicate, chromosome.
select.downsample.ids <- function(
    data, downsample.size, sample.id.column, grouping.columns, seed
  ) {
  candidates <- data %>%
    filter(data.type != "Empirical") %>%
    distinct(across(all_of(c(grouping.columns, sample.id.column))))
  sizes <- candidates %>% count(across(all_of(grouping.columns)), name = "n")
  if (any(sizes$n < downsample.size)) {
    stop("A simulation group contains fewer than ", downsample.size,
      " candidates")
    }
  set.seed(seed)
  selected.ids <- candidates %>%
    group_by(across(all_of(grouping.columns))) %>%
    slice_sample(n = downsample.size, replace = FALSE) %>%
    ungroup() %>%
    arrange(across(all_of(grouping.columns)), .data[[sample.id.column]])
  return(selected.ids)
  }


# add complete and selected simulation sample sets without duplicating
# empirical rows.
apply.downsample.ids <- function(
    data, selected.ids, sample.id.column, grouping.columns
  ) {
  duplicates <- selected.ids %>%
    count(across(all_of(c(grouping.columns, sample.id.column)))) %>%
    filter(n != 1L)
  if (nrow(duplicates)) stop("Selected sample IDs are not unique")
  full <- mutate(data, sample.set = "full")
  downsampled <- data %>%
    filter(data.type != "Empirical") %>%
    inner_join(selected.ids, by = c(grouping.columns, sample.id.column)) %>%
    mutate(sample.set = "downsampled")
  return(bind_rows(full, downsampled))
  }


# evaluate ancestry statistics and reject non-finite values.
calculate.ancestry.statistics <- function(values) {
  statistics <- c(
    mean = mean(values),
    sd = sd(values)
    )
  if (any(!is.finite(statistics))) {
    stop("Ancestry mean and SD must be finite")
    }
  return(statistics)
  }


# create mean and sd inputs for bootstrap comparison tables.
summarize.ancestry.comparison <- function(data) {
  summary <- data %>%
    filter(!is.na(afr.q)) %>%
    group_by(rep, chrom, data.type, method, sample.set) %>%
    group_modify(function(group, key) {
      return(as_tibble_row(calculate.ancestry.statistics(group$afr.q)))
      }) %>%
    ungroup()
  return(summary)
  }


# summarize complete simulation replicates with percentile intervals.
summarize.bootstrap.interval <- function(data, grouping.columns, value.column) {
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


# retain one active ancestry method for each simulation source.
filter.bootstrap.ancestry.simulation <- function(data) {
  simulation <- data %>%
    filter(
      data.type != "Empirical", sample.set == "full",
      (data.type %in% ANCESTRY.SOURCE.LEVELS[c(1, 3)] & method == "tspop") |
        (data.type %in% ANCESTRY.SOURCE.LEVELS[c(2, 4)] &
          method == ANCESTRY.PLOT.EMPIRICAL.METHOD)
      )
  return(simulation)
  }


# select empirical-sized adx samples independently per simulation group.
select.bootstrap.ancestry.ids <- function(
    data, downsample.size, sample.id.column, seed
  ) {
  candidates <- filter.bootstrap.ancestry.simulation(data) %>%
    select(data.type, rep, chrom, all_of(sample.id.column), afr.q)
  if (any(!is.finite(candidates$afr.q))) {
    stop("Simulation ancestry values must be finite")
    }
  duplicates <- candidates %>%
    count(data.type, rep, chrom, .data[[sample.id.column]]) %>%
    filter(n != 1L)
  if (nrow(duplicates)) stop("Simulation ancestry IDs must be unique")
  sizes <- candidates %>% count(data.type, rep, chrom, name = "available")
  if (any(sizes$available < downsample.size)) {
    stop("A simulation ancestry group contains fewer than ", downsample.size,
      " unique IDs")
    }
  set.seed(seed)
  selected <- candidates %>%
    group_by(data.type, rep, chrom) %>%
    slice_sample(n = downsample.size, replace = FALSE) %>%
    ungroup() %>%
    select(data.type, rep, chrom, all_of(sample.id.column))
  return(selected)
  }


# build simulation intervals and empirical bootstrap statistic summaries.
summarize.bootstrap.ancestry <- function(
    data, downsample.size, seed, replicates
  ) {
  selected <- select.bootstrap.ancestry.ids(
    data, downsample.size, "sample_id", seed
    )
  simulation <- filter.bootstrap.ancestry.simulation(data) %>%
    inner_join(selected, by = c("data.type", "rep", "chrom", "sample_id")) %>%
    group_by(data.type, rep, chrom) %>%
    group_modify(function(group, key) {
      return(as_tibble_row(calculate.ancestry.statistics(group$afr.q)))
      }) %>%
    ungroup() %>%
    pivot_longer(c(mean, sd), names_to = "stat", values_to = "value") %>%
    summarize.bootstrap.interval(c("data.type", "rep", "chrom", "stat"),
      "value")
  empirical <- data %>%
    filter(
      data.type == "Empirical", sample.set == "full",
      method == ANCESTRY.PLOT.EMPIRICAL.METHOD
      )
  if (any(!is.finite(empirical$afr.q))) {
    stop("Empirical ancestry values must be finite")
    }
  set.seed(seed)
  empirical <- empirical %>%
    group_by(data.type, chrom) %>%
    group_modify(function(group, key) {
      observed <- calculate.ancestry.statistics(group$afr.q)
      estimates <- map_dfr(names(observed), function(statistic) {
        values <- replicate(replicates, {
          calculate.ancestry.statistics(sample(
            group$afr.q, nrow(group), replace = TRUE
            ))[[statistic]]
          })
        tibble(
          stat = statistic, mean = observed[[statistic]],
          lower = quantile(values, 0.025, names = FALSE),
          upper = quantile(values, 0.975, names = FALSE),
          replicate.count = as.integer(replicates)
          )
        })
      return(estimates)
      }) %>%
    ungroup()
  return(bind_rows(simulation, empirical))
  }


# summarize simulation and genome-wide asw bootstrap histogram fractions.
summarize.bootstrap.histograms <- function(data, breaks, seed, replicates) {
  selected <- select.bootstrap.ancestry.ids(data, ANCESTRY.DOWNSAMPLE.SIZE,
    "sample_id", seed)
  simulation <- filter.bootstrap.ancestry.simulation(data) %>%
    inner_join(selected, by = c("data.type", "rep", "chrom", "sample_id")) %>%
    group_by(data.type, rep, chrom) %>%
    group_modify(function(group, key) {
      counts <- hist(group$afr.q, breaks = breaks, plot = FALSE)$counts
      return(tibble(bin = seq_along(counts), fraction = counts / sum(counts)))
      }) %>%
    ungroup() %>%
    summarize.bootstrap.interval(c("data.type", "rep", "chrom", "bin"),
      "fraction")
  empirical <- data %>%
    filter(
      data.type == "Empirical", sample.set == "full", chrom == "all",
      method == ANCESTRY.PLOT.EMPIRICAL.METHOD
      ) %>%
    group_by(data.type, chrom) %>%
    group_modify(function(group, key) {
      set.seed(seed)
      draws <- replicate(replicates, {
        counts <- hist(sample(group$afr.q, nrow(group), replace = TRUE),
          breaks = breaks, plot = FALSE)$counts
        counts / sum(counts)
        })
      return(tibble(
        bin = seq_len(nrow(draws)), mean = rowMeans(draws),
        lower = apply(draws, 1, quantile, 0.025),
        upper = apply(draws, 1, quantile, 0.975),
        replicate.count = replicates
        ))
      }) %>%
    ungroup()
  bins <- tibble(
    bin = seq_len(length(breaks) - 1L),
    xmid = head(breaks, -1L) + diff(breaks) / 2
    )
  return(bind_rows(simulation, empirical) %>% left_join(bins, by = "bin"))
  }


# construct selected-chromosome bars with a genome-wide asw reference.
make.bootstrap.ancestry.bar.plot <- function(data, data.types, title) {
  plotted <- data %>%
    filter(
      as.character(data.type) %in% data.types,
      chrom %in% ANCESTRY.SELECTED.CHROMOSOMES
      ) %>%
    mutate(data.type = factor(
      as.character(data.type),
      levels = order.active.levels(data.type, ANCESTRY.SOURCE.LEVELS)
      ))
  genome <- data %>%
    filter(data.type == "Empirical", chrom == "all", stat %in% plotted$stat)
  dodge <- position_dodge(width = ANCESTRY.CATEGORICAL.BAR.DODGE)
  plot <- ggplot(plotted, aes(chrom, mean, fill = data.type)) +
    geom_rect(
      data = genome, aes(ymin = lower, ymax = upper),
      xmin = -Inf, xmax = Inf, inherit.aes = FALSE,
      fill = ANCESTRY.PLOT.STYLES$empirical.colors[[ANCESTRY.PLOT.EMPIRICAL.METHOD]],
      alpha = 0.15
      ) +
    geom_col(
      position = dodge, width = ANCESTRY.CATEGORICAL.BAR.WIDTH,
      linewidth = ANCESTRY.CATEGORICAL.BAR.LINEWIDTH, color = "black"
      ) +
    geom_errorbar(
      aes(ymin = lower, ymax = upper), position = dodge, width = 0,
      linewidth = ANCESTRY.CATEGORICAL.BAR.LINEWIDTH, na.rm = TRUE
      ) +
    geom_hline(
      data = genome, aes(yintercept = mean), linetype = "longdash",
      color = ANCESTRY.PLOT.STYLES$empirical.colors[[ANCESTRY.PLOT.EMPIRICAL.METHOD]],
      linewidth = 1
      ) +
    facet_wrap(
      vars(stat), scales = "free_y", nrow = 1,
      labeller = labeller(stat = c(mean = "Mean", sd = "SD"))
      ) +
    facetted_pos_scales(
      y = list(
        stat == "mean" ~ scale_y_continuous(limits = c(0.70, 0.95)),
        stat == "sd" ~ scale_y_continuous(limits = c(0, 0.25))
        )
      ) +
    scale_fill_manual(values = ANCESTRY.PLOT.STYLES$colors,
      labels = ANCESTRY.PLOT.STYLES$labels) +
    labs(
      title = title, subtitle = ancestry.inference.subtitle(), x = "Chromosome",
      y = "African ancestry", fill = NULL
      )
  return(apply.standard.plot.theme(plot))
  }


# construct one simulated-chromosome and genome-wide asw histogram.
make.bootstrap.ancestry.histogram.plot <- function(
    data, data.types, simulated.chromosome, title
  ) {
  if (length(simulated.chromosome) != 1L ||
      !is.finite(as.numeric(simulated.chromosome))) {
    stop("Ancestry histogram requires one finite simulated chromosome")
    }
  simulated.chromosome <- as.character(simulated.chromosome)
  plotted <- data %>%
    filter(
      as.character(data.type) %in% data.types,
      (as.character(chrom) == "all" & data.type == "Empirical") |
        (as.character(chrom) == simulated.chromosome &
          data.type != "Empirical")
      ) %>%
    mutate(data.type = factor(
      as.character(data.type),
      levels = order.active.levels(data.type, ANCESTRY.SOURCE.LEVELS)
      ))
  dodge <- position_dodge(width = diff(ANCESTRY.HISTOGRAM.BREAKS)[1])
  plot <- ggplot(plotted, aes(xmid, mean, fill = data.type)) +
    geom_col(
      position = dodge,
      width = diff(ANCESTRY.HISTOGRAM.BREAKS)[1] * ANCESTRY.DENSE.BAR.WIDTH.MULTIPLIER,
      color = "black", linewidth = ANCESTRY.DENSE.BAR.LINEWIDTH
      ) +
    geom_errorbar(
      aes(ymin = lower, ymax = upper), position = dodge, width = 0,
      linewidth = ANCESTRY.DENSE.BAR.LINEWIDTH, na.rm = TRUE
      ) +
    scale_fill_manual(values = ANCESTRY.PLOT.STYLES$colors,
      labels = ANCESTRY.PLOT.STYLES$labels) +
    labs(
      title = title, x = "African ancestry", y = "Fraction of individuals",
      fill = NULL, subtitle = ancestry.inference.subtitle()
      )
  return(apply.standard.plot.theme(plot))
  }


# validate summary and asw values used by bootstrap difference tables.
prepare.ancestry.bootstrap.comparison.data <- function(
    summary.data, individual.data, empirical.method, sample.set.input,
    chromosomes
  ) {
  chromosomes <- as.character(chromosomes)
  if (!setequal(chromosomes, as.character(1:22)) ||
      length(unique(chromosomes)) != 22L) {
    stop("Bootstrap ancestry comparisons require chromosomes 1-22")
    }
  required.summary <- c(
    "data.type", "method", "sample.set", "chrom", "rep", "mean", "sd"
    )
  missing.summary <- setdiff(required.summary, names(summary.data))
  if (length(missing.summary)) {
    stop("Bootstrap ancestry summaries are missing columns: ",
      paste(missing.summary, collapse = ", "))
    }
  source.config <- tribble(
    ~source, ~data.type, ~method,
    "TC", "Simulation_2T12Consistent", "tspop",
    "TCD", "Simulation_2T12Consistent_simDown", empirical.method,
    "LG", "Simulation_largeGrowth", "tspop",
    "LGD", "Simulation_largeGrowth_simDown", empirical.method
    )
  simulation <- summary.data %>%
    mutate(chrom = as.character(chrom)) %>%
    inner_join(source.config, by = c("data.type", "method")) %>%
    filter(sample.set == sample.set.input, chrom %in% chromosomes) %>%
    select(source, chrom, rep, mean, sd)
  if (any(!is.finite(unlist(simulation[c("mean", "sd")]))) ||
      any(!is.finite(simulation$rep))) {
    stop("Bootstrap ancestry comparison summaries must be finite")
    }
  duplicates <- simulation %>%
    count(source, chrom, rep, name = "rows") %>%
    filter(rows != 1L)
  if (nrow(duplicates)) {
    first.duplicate <- duplicates[1, ]
    stop(
      "Bootstrap ancestry comparison summaries contain duplicated replicate ",
      "IDs for ", first.duplicate$source, " chromosome ",
      first.duplicate$chrom
      )
    }
  counts <- simulation %>%
    count(source, chrom, name = "replicate.count") %>%
    complete(
      source = source.config$source, chrom = chromosomes,
      fill = list(replicate.count = 0L)
      ) %>%
    filter(replicate.count != 50L)
  if (nrow(counts)) {
    first.count <- counts[1, ]
    stop(
      "Each bootstrap ancestry comparison source and chromosome must contain ",
      "exactly 50 replicate IDs; found ", first.count$replicate.count,
      " for ", first.count$source, " chromosome ", first.count$chrom
      )
    }
  pwalk(
    filter(ANCESTRY.BOOTSTRAP.CONTRASTS, paired),
    function(contrast, left.source, right.source, paired) {
      for (chromosome in chromosomes) {
        left.ids <- simulation %>%
          filter(source == left.source, chrom == chromosome) %>% pull(rep)
        right.ids <- simulation %>%
          filter(source == right.source, chrom == chromosome) %>% pull(rep)
        if (!setequal(left.ids, right.ids)) {
          stop(contrast, " chromosome ", chromosome,
            " requires matching replicate IDs", call. = FALSE)
          }
        }
      }
    )
  required.individual <- c(
    "data.type", "role", "method", "chrom", "sample_id", "afr.q"
    )
  missing.individual <- setdiff(required.individual, names(individual.data))
  if (length(missing.individual)) {
    stop("Bootstrap empirical ancestry data are missing columns: ",
      paste(missing.individual, collapse = ", "))
    }
  empirical <- individual.data %>%
    mutate(
      chrom = as.character(chrom),
      sample_id = as.character(vcf_sample_id) # copy vcf_sample_id into sample_id
      ) %>%
    filter(
      data.type == "Empirical", role == "ASW", method == empirical.method,
      chrom %in% c(chromosomes, "all")
      )
  if ("sample.set" %in% names(empirical)) {
    empirical <- filter(empirical, sample.set == "full")
    }
  if (any(!is.finite(empirical$afr.q))) {
    stop("Bootstrap empirical ASW ancestry values must be finite")
    }
  if (any(is.na(empirical$sample_id) |
      !nzchar(trimws(empirical$sample_id)))) {
    stop("Bootstrap empirical ASW ancestry requires a valid sample_id")
    }
  empirical.counts <- empirical %>%
    count(chrom, name = "individual.count") %>%
    complete(
      chrom = c(chromosomes, "all"), fill = list(individual.count = 0L)
      ) %>%
    filter(individual.count < 2L)
  if (nrow(empirical.counts)) {
    first.count <- empirical.counts[1, ]
    stop("Bootstrap empirical ASW ancestry requires at least two individuals ",
      "for chromosome ", first.count$chrom)
    }
  duplicates <- empirical %>%
    count(chrom, sample_id, name = "rows") %>%
    filter(rows != 1L)
  if (nrow(duplicates)) {
    first.duplicate <- duplicates[1, ]
    stop(
      "Bootstrap empirical ASW sample IDs contain duplicated chromosome/ID: ",
      first.duplicate$chrom, " / ", first.duplicate$sample_id
      )
    }
  return(list(simulation = simulation, empirical = empirical))
  }


# summarize a bootstrap difference and nominal and corrected intervals.
summarize.ancestry.bootstrap.difference <- function(
    draws, difference, family.size
  ) {
  if (any(!is.finite(draws)) || !is.finite(difference)) {
    stop("Bootstrap ancestry comparison differences must be finite")
    }
  nominal <- quantile(draws, c(0.025, 0.975), names = FALSE)
  corrected <- quantile(
    draws,
    c(0.05 / (2 * family.size), 1 - 0.05 / (2 * family.size)),
    names = FALSE
    )
  return(tibble(
    difference = difference,
    ci.95.lower = nominal[[1]], ci.95.upper = nominal[[2]],
    bonferroni.ci.lower = corrected[[1]],
    bonferroni.ci.upper = corrected[[2]]
    ))
  }


# resample replicate-level simulation values for one contrast.
bootstrap.ancestry.simulation.difference <- function(
    left.values, right.values, paired, bootstrap.replicates
  ) {
  if (paired) {
    draws <- replicate(bootstrap.replicates, {
      indices <- sample(seq_along(left.values), length(left.values),
        replace = TRUE)
      mean(left.values[indices] - right.values[indices])
      })
    } else {
    draws <- replicate(bootstrap.replicates, {
      mean(sample(left.values, length(left.values), replace = TRUE)) -
        mean(sample(right.values, length(right.values), replace = TRUE))
      })
    }
  return(list(
    draws = draws,
    difference = mean(left.values) - mean(right.values)
    ))
  }


# resample a simulation summary and asw individuals for one contrast.
bootstrap.ancestry.empirical.difference <- function(
    simulation.values, empirical.values, statistic, bootstrap.replicates
  ) {
  draws <- replicate(bootstrap.replicates, {
    mean(sample(simulation.values, length(simulation.values), replace = TRUE)) -
      calculate.ancestry.statistics(sample(
        empirical.values, length(empirical.values), replace = TRUE
        ))[[statistic]]
    })
  difference <- mean(simulation.values) -
    calculate.ancestry.statistics(empirical.values)[[statistic]]
  return(list(draws = draws, difference = difference))
  }


# build one table against chromosome or genome asw values.
make.ancestry.bootstrap.comparison.table <- function(
    simulation, empirical, contrasts, chromosomes, statistic, reference.chrom,
    reference.scope, family.size, bootstrap.replicates, seed = ANCESTRY.RANDOM.SEED
  ) {
  work <- crossing(
    chromosome = chromosomes, contrast.index = seq_len(nrow(contrasts))
    ) %>%
    mutate(
      contrast.data = map(contrast.index, ~ contrasts[.x, ]),
      reference.chromosome = map_chr(chromosome, reference.chrom)
      ) %>%
    mutate(
      left.values = pmap(list(chromosome, contrast.data),
        function(chromosome, contrast.data) simulation %>% filter(
          source == contrast.data$left.source, chrom == .env$chromosome
          ) %>% arrange(rep) %>% pull(all_of(statistic))),
      right.values = pmap(list(chromosome, contrast.data),
        function(chromosome, contrast.data) simulation %>% filter(
          source == contrast.data$right.source, chrom == .env$chromosome
          ) %>% arrange(rep) %>% pull(all_of(statistic))),
      empirical.values = map(reference.chromosome,
        ~ empirical %>% filter(chrom == .x) %>% pull(afr.q))
      )
  table <- bind_rows(bootstrap.parallel.pmap(
    work,
    function(chromosome, contrast.index, contrast.data, reference.chromosome,
             left.values, right.values, empirical.values) {
      bootstrap <- if (contrast.data$right.source == "Emp") {
        bootstrap.ancestry.empirical.difference(
          left.values, empirical.values, statistic, bootstrap.replicates
          )
        } else {
        bootstrap.ancestry.simulation.difference(
          left.values, right.values, contrast.data$paired, bootstrap.replicates
          )
        }
      bind_cols(
        tibble(statistic = statistic, contrast = contrast.data$contrast,
          chromosome = chromosome, reference.scope = reference.scope),
        summarize.ancestry.bootstrap.difference(
          bootstrap$draws, bootstrap$difference, family.size
          )
        )
      },
    seed = seed
    ))
  return(table)
  }


# make deterministic chromosome and genome-wide asw comparison tables.
make.ancestry.bootstrap.comparison.tables <- function(
    summary.data, individual.data, empirical.method = ANCESTRY.PLOT.EMPIRICAL.METHOD,
    sample.set.input = ANCESTRY.PLOT.SAMPLE.SET, chromosomes = ANCESTRY.CHROMOSOMES,
    bootstrap.replicates = ANCESTRY.BOOTSTRAP.REPLICATES, seed = ANCESTRY.RANDOM.SEED
  ) {
  if (!is.numeric(bootstrap.replicates) || length(bootstrap.replicates) != 1L ||
      !is.finite(bootstrap.replicates) || bootstrap.replicates < 1L ||
      bootstrap.replicates %% 1L != 0) {
    stop("Bootstrap comparison replicate count must be a positive integer")
    }
  data <- prepare.ancestry.bootstrap.comparison.data(
    summary.data, individual.data, empirical.method, sample.set.input,
    chromosomes
    )
  chromosomes <- as.character(chromosomes)
  genome.contrasts <- filter(
    ANCESTRY.BOOTSTRAP.CONTRASTS, right.source == "Emp"
    )
  chromosome.tables <- map2_dfr(c("mean", "sd"), c(0L, 1L),
      function(statistic, seed.offset) {
    make.ancestry.bootstrap.comparison.table(
      data$simulation, data$empirical, ANCESTRY.BOOTSTRAP.CONTRASTS,
      chromosomes, statistic, identity, "chromosome",
      ANCESTRY.BOOTSTRAP.CHROMOSOME.FAMILY.SIZE, bootstrap.replicates,
      seed + seed.offset
      )
    })
  genome.tables <- map2_dfr(c("mean", "sd"), c(2L, 3L),
      function(statistic, seed.offset) {
    make.ancestry.bootstrap.comparison.table(
      data$simulation, data$empirical, genome.contrasts, chromosomes,
      statistic, function(chromosome) "all", "genome",
      ANCESTRY.BOOTSTRAP.GENOME.ASW.FAMILY.SIZE, bootstrap.replicates,
      seed + seed.offset
      )
    })
  if (nrow(chromosome.tables) != 308L) {
    stop(
      "Bootstrap chromosome ancestry comparison table row count is incorrect"
      )
    }
  if (nrow(genome.tables) != 176L) {
    stop(
      "Bootstrap genome-ASW ancestry comparison table row count is incorrect"
      )
    }
  return(list(
    chromosome.comparisons = chromosome.tables,
    chromosome.vs.genome.asw = genome.tables
    ))
  }


# persist the two deterministic bootstrap comparison tables.
write.ancestry.bootstrap.comparison.tables <- function(
    tables, output.directory
  ) {
  dir.create(output.directory, recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(
    tables$chromosome.comparisons,
    file.path(
      output.directory,
      ancestry.output.filename(
        "ancestry.statistics.bootstrap.chromosome.comparisons", "csv"
        )
      )
    )
  readr::write_csv(
    tables$chromosome.vs.genome.asw,
    file.path(
      output.directory,
      ancestry.output.filename(
        "ancestry.statistics.bootstrap.chromosome_vs_genome_asw", "csv"
        )
      )
    )
  }


# prepare one contrast family for nominal or Bonferroni interval plotting.
prepare.ancestry.bootstrap.contrast.plot.data <- function(
    tables, contrast.family, interval.type = c("95", "bonferroni"),
    reference.scope = c("chromosome", "genome"),
    statistics = c("mean", "sd")
  ) {
  interval.type <- match.arg(interval.type)
  if (contrast.family == "empirical") {
    reference.scope <- match.arg(reference.scope)
    }
  contrast.config <- list(
    empirical = list(
      contrasts = c("TC-Emp", "TCD-Emp", "LG-Emp", "LGD-Emp"),
      chromosome.levels = ANCESTRY.CHROMOSOMES
      ),
    simulation = list(
      contrasts = c("TC-TCD", "LG-LGD", "TC-LG"),
      chromosome.levels = ANCESTRY.CHROMOSOMES
      )
    )
  if (!contrast.family %in% names(contrast.config)) {
    stop("Unsupported bootstrap contrast family: ", contrast.family)
    }
  config <- contrast.config[[contrast.family]]
  if (!all(statistics %in% c("mean", "sd"))) {
    stop("Unsupported ancestry bootstrap statistic")
    }
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
  data <- if (contrast.family == "empirical") {
    if (reference.scope == "chromosome") {
      tables$chromosome.comparisons
      } else {
      tables$chromosome.vs.genome.asw
      }
    } else {
    tables$chromosome.comparisons
    }
  plotted <- data %>%
    filter(contrast %in% config$contrasts, statistic %in% statistics) %>%
    mutate(
      ci.lower = .data[[lower.column]],
      ci.upper = .data[[upper.column]],
      chromosome = factor(chromosome, levels = config$chromosome.levels),
      contrast = factor(contrast, levels = config$contrasts),
      statistic = factor(statistic, levels = statistics),
      significant = ci.lower > 0 | ci.upper < 0,
      point.color = if_else(significant, "red", "black"),
      contrast.color = unname(ANCESTRY.PLOT.STYLES$contrast.colors[as.character(
        contrast
        )]),
      facet.label = factor(NA_character_)
      ) %>%
    group_by(chromosome) %>%
    mutate(plot.x = as.numeric(chromosome)) %>%
    ungroup() %>%
    arrange(statistic, chromosome, contrast)
  return(plotted)
  }


# construct a faceted bootstrap contrast plot from prepared comparison data.
make.ancestry.bootstrap.contrast.plot <- function(
    data, title, facet.statistics = TRUE
  ) {
  contrasts <- levels(data$contrast)
  active.contrasts <- contrasts[
    contrasts %in% unique(as.character(data$contrast))
    ]
  plot <- ggplot(
    data,
    aes(plot.x, difference, group = contrast)
    ) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1) +
    geom_errorbar(
      aes(ymin = ci.lower, ymax = ci.upper), width = 0,
      color = "black", linewidth = ANCESTRY.CATEGORICAL.BAR.LINEWIDTH,
      show.legend = FALSE
      ) +
    geom_point(
      aes(fill = contrast, color = point.color),
      shape = 21, size = 3, stroke = ANCESTRY.CATEGORICAL.BAR.LINEWIDTH
      ) +
    scale_color_manual(
      values = c(red = "red", black = "black", ANCESTRY.PLOT.STYLES$contrast.colors),
      guide = "none") +
    scale_fill_manual(
      values = ANCESTRY.PLOT.STYLES$contrast.colors, breaks = active.contrasts,
      labels = ANCESTRY.PLOT.STYLES$contrast.labels[active.contrasts]
      ) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
    labs(
      x = "Chromosome", y = "Difference", color = NULL, fill = NULL,
      title = title,
      subtitle = ancestry.inference.subtitle()
      ) +
    guides(
      fill = guide_legend(
        override.aes = list(
          shape = 21,
          fill = unname(ANCESTRY.PLOT.STYLES$contrast.colors[active.contrasts]),
          color = "black"
          )
        )
      ) +
    scale_x_continuous(
      breaks = seq_along(levels(data$chromosome)),
      labels = c("1", rep("", 9), "11", rep("", 10), "22")
      ) +
    theme(legend.title = element_blank())
    
  if (facet.statistics) {
    if (all(is.na(data$facet.label))) {
      plot <- plot + facet_grid(
        cols = vars(contrast), rows = vars(statistic), scales = "free_y",
        labeller = labeller(
          contrast = ANCESTRY.PLOT.STYLES$contrast.labels,
          statistic = c(mean = "Mean", sd = "SD"))
        )
      } else {
      plot <- plot + facet_wrap(
        cols = vars(contrast), rowsvars(facet.label), ncol = 1, 
        scales = "free_y",
        strip.position = "right"
        )
      }
    }
  return(apply.standard.plot.theme(plot))
  }


# analysis ----


# read tc/lg truth, configured tcd/lgd inference, and empirical inference.
ancestry.sim.tc.tspop.data <- read.ancestry.family(
  ANCESTRY.SIM.TC.DATA.DIR, "ancestry.chr{chrom}.parquet", ANCESTRY.CHROMOSOMES,
  "Simulation_2T12Consistent", "tspop", 0L, "tspop"
  )
ancestry.simDown.tc.inference.data <- read.ancestry.family(
  ANCESTRY.SIMDOWN.TC.DATA.DIR, ancestry.inference.file.family(ANCESTRY.PLOT.EMPIRICAL.METHOD),
  ANCESTRY.CHROMOSOMES, "Simulation_2T12Consistent_simDown", ANCESTRY.PLOT.EMPIRICAL.METHOD,
  ANCESTRY.SIMULATION.K, ANCESTRY.PLOT.EMPIRICAL.METHOD
  )
ancestry.sim.lg.tspop.data <- read.ancestry.family(
  ANCESTRY.SIM.LG.DATA.DIR, "ancestry.chr{chrom}.parquet", ANCESTRY.CHROMOSOMES,
  "Simulation_largeGrowth", "tspop", 0L, "tspop"
  )
ancestry.simDown.lg.inference.data <- read.ancestry.family(
  ANCESTRY.SIMDOWN.LG.DATA.DIR, ancestry.inference.file.family(ANCESTRY.PLOT.EMPIRICAL.METHOD),
  ANCESTRY.CHROMOSOMES, "Simulation_largeGrowth_simDown", ANCESTRY.PLOT.EMPIRICAL.METHOD,
  ANCESTRY.SIMULATION.K, ANCESTRY.PLOT.EMPIRICAL.METHOD
  )
ancestry.empirical.file.family <- ancestry.inference.file.family(
  ANCESTRY.PLOT.EMPIRICAL.METHOD
  )
ancestry.empirical.genome.file <- str_replace(
  ancestry.empirical.file.family, "\\.chr\\{chrom\\}", ""
  )
ancestry.empirical.inference.data <- read.ancestry.family(
  ANCESTRY.EMPIRICAL.DATA.DIR, ancestry.empirical.file.family, ANCESTRY.CHROMOSOMES, "Empirical",
  ANCESTRY.PLOT.EMPIRICAL.METHOD, ANCESTRY.EMPIRICAL.K, "Empirical", include.genome = TRUE,
  genome.file.family = ancestry.empirical.genome.file
  )

# orient ancestry components and add reproducible selected simulation rows.
ancestry.individual.data <- bind_rows(
  ancestry.sim.tc.tspop.data, ancestry.simDown.tc.inference.data,
  ancestry.sim.lg.tspop.data, ancestry.simDown.lg.inference.data,
  ancestry.empirical.inference.data
  ) %>%
  apply.ancestry.source.contract() %>%
  orient.ancestry.components(c("rep", "chrom", "data.type", "method"))
ancestry.downsample.ids <- select.downsample.ids(
  ancestry.individual.data, ANCESTRY.DOWNSAMPLE.SIZE, "sample_id",
  c("data.type", "rep", "chrom"), ANCESTRY.RANDOM.SEED
  )
ancestry.individual.data <- apply.downsample.ids(
  ancestry.individual.data, ancestry.downsample.ids, "sample_id",
  c("data.type", "rep", "chrom")
  )

# build and write the two active bootstrap comparison tables.
ancestry.summary.data <- summarize.ancestry.comparison(ancestry.individual.data)
ancestry.bootstrap.comparison.tables <-
  make.ancestry.bootstrap.comparison.tables(
  ancestry.summary.data, ancestry.individual.data,
  empirical.method = ANCESTRY.PLOT.EMPIRICAL.METHOD,
  sample.set.input = ANCESTRY.PLOT.SAMPLE.SET, chromosomes = ANCESTRY.CHROMOSOMES,
  bootstrap.replicates = ANCESTRY.BOOTSTRAP.REPLICATES, seed = ANCESTRY.RANDOM.SEED
  )
write.ancestry.bootstrap.comparison.tables(
  ancestry.bootstrap.comparison.tables, ANCESTRY.OUTPUT.DIR
  )


# construct the four active bottom bootstrap figures.
ancestry.bootstrap.summary <- summarize.bootstrap.ancestry(
  ancestry.individual.data, ANCESTRY.DOWNSAMPLE.SIZE, ANCESTRY.RANDOM.SEED,
  ANCESTRY.BOOTSTRAP.REPLICATES
  )
ancestry.bootstrap.histograms <- summarize.bootstrap.histograms(
  ancestry.individual.data, ANCESTRY.HISTOGRAM.BREAKS, ANCESTRY.RANDOM.SEED,
  ANCESTRY.BOOTSTRAP.REPLICATES
  )


# plotting ----
ancestry.bootstrap.tcd.1kg.bar <- make.bootstrap.ancestry.bar.plot(
  ancestry.bootstrap.summary,
  c("Simulation_2T12Consistent_simDown", "Empirical"),
  "African ancestry: TCD and ASW"
)
ancestry.bootstrap.all.datatypes.adx.asw.bar <-
  make.bootstrap.ancestry.bar.plot(
    ancestry.bootstrap.summary, ANCESTRY.SOURCE.LEVELS,
    "African ancestry: all ADX sources and ASW"
  )
ancestry.bootstrap.histogram.chromosomes <- as.character(1:22)
ancestry.bootstrap.tcd.1kg.histograms <- setNames(
  lapply(ancestry.bootstrap.histogram.chromosomes, function(chromosome) {
    make.bootstrap.ancestry.histogram.plot(
      ancestry.bootstrap.histograms,
      c("Simulation_2T12Consistent_simDown", "Empirical"), chromosome,
      paste0(
        "Chromosome ", chromosome,
        " simulations and genome-wide ASW: TCD and ASW"
      )
    )
  }),
  paste0("chr", ancestry.bootstrap.histogram.chromosomes)
)
ancestry.bootstrap.all.datatypes.adx.asw.histograms <- setNames(
  lapply(ancestry.bootstrap.histogram.chromosomes, function(chromosome) {
    make.bootstrap.ancestry.histogram.plot(
      ancestry.bootstrap.histograms, ANCESTRY.SOURCE.LEVELS, chromosome,
      paste0(
        "Chromosome ", chromosome,
        " simulations and genome-wide ASW: all ADX sources and ASW"
      )
    )
  }),
  paste0("chr", ancestry.bootstrap.histogram.chromosomes)
)

# retain chromosome-1 aliases used by existing downstream consumers.
ancestry.bootstrap.tcd.1kg.histogram <-
  ancestry.bootstrap.tcd.1kg.histograms$chr1
ancestry.bootstrap.all.datatypes.adx.asw.histogram <-
  ancestry.bootstrap.all.datatypes.adx.asw.histograms$chr1
ancestry.bootstrap.empirical.95.chromosome.comparisons <-
  make.ancestry.bootstrap.contrast.plot(
    prepare.ancestry.bootstrap.contrast.plot.data(
      ancestry.bootstrap.comparison.tables, "empirical", "95", "chromosome"
    ),
    "Bootstrap ancestry differences: chromosome by chromosome"
  )
ancestry.bootstrap.empirical.95.genome.comparisons <-
  make.ancestry.bootstrap.contrast.plot(
    prepare.ancestry.bootstrap.contrast.plot.data(
      ancestry.bootstrap.comparison.tables, "empirical", "95", "genome"
    ),
    "Bootstrap ancestry differences: chromosome by whole-genome ASW"
  )
ancestry.bootstrap.empirical.bonferroni.chromosome.comparisons <-
  make.ancestry.bootstrap.contrast.plot(
    prepare.ancestry.bootstrap.contrast.plot.data(
      ancestry.bootstrap.comparison.tables, "empirical", "bonferroni",
      "chromosome"
    ),
    "Bonferroni bootstrap ancestry differences: chromosome by chromosome"
  )
ancestry.bootstrap.empirical.bonferroni.genome.comparisons <-
  make.ancestry.bootstrap.contrast.plot(
    prepare.ancestry.bootstrap.contrast.plot.data(
      ancestry.bootstrap.comparison.tables, "empirical", "bonferroni",
      "genome"
    ),
    "Bonferroni bootstrap ancestry differences: chromosome by whole-genome ASW"
  )
ancestry.bootstrap.simulation.95.comparisons <-
  make.ancestry.bootstrap.contrast.plot(
    prepare.ancestry.bootstrap.contrast.plot.data(
      ancestry.bootstrap.comparison.tables, "simulation", "95"
    ),
    "Bootstrap ancestry differences: simulation contrasts"
  )
ancestry.bootstrap.simulation.bonferroni.comparisons <-
  make.ancestry.bootstrap.contrast.plot(
    prepare.ancestry.bootstrap.contrast.plot.data(
      ancestry.bootstrap.comparison.tables, "simulation", "bonferroni"
    ),
    "Bonferroni bootstrap ancestry differences: simulation contrasts"
  )

# save each active figure before explicit printing at the script end.
dir.create(ANCESTRY.OUTPUT.DIR, recursive = TRUE, showWarnings = FALSE)
saveRDS(ancestry.bootstrap.tcd.1kg.bar, file.path(
  ANCESTRY.OUTPUT.DIR, ancestry.output.filename("ancestry.bootstrap.tcd.1kg.bar", "rds")
  ))
saveRDS(ancestry.bootstrap.all.datatypes.adx.asw.bar, file.path(
  ANCESTRY.OUTPUT.DIR, ancestry.output.filename(
    "ancestry.bootstrap.all.datatypes.adx.asw.bar", "rds"
    )
  ))
saveRDS(ancestry.bootstrap.tcd.1kg.histogram, file.path(
  ANCESTRY.OUTPUT.DIR, ancestry.output.filename(
    "ancestry.bootstrap.tcd.1kg.histogram", "rds"
    )
  ))
saveRDS(ancestry.bootstrap.all.datatypes.adx.asw.histogram, file.path(
  ANCESTRY.OUTPUT.DIR, ancestry.output.filename(
    "ancestry.bootstrap.all.datatypes.adx.asw.histogram", "rds"
    )
  ))
saveRDS(ancestry.bootstrap.tcd.1kg.histograms, file.path(
  ANCESTRY.OUTPUT.DIR, ancestry.output.filename(
    "ancestry.bootstrap.tcd.1kg.histograms", "rds"
    )
  ))
saveRDS(ancestry.bootstrap.all.datatypes.adx.asw.histograms, file.path(
  ANCESTRY.OUTPUT.DIR, ancestry.output.filename(
    "ancestry.bootstrap.all.datatypes.adx.asw.histograms", "rds"
    )
  ))
saveRDS(ancestry.bootstrap.empirical.95.chromosome.comparisons, file.path(
  ANCESTRY.OUTPUT.DIR, ancestry.output.filename(
    "ancestry.bootstrap.empirical.95.chromosome.comparisons", "rds"
    )
  ))
saveRDS(ancestry.bootstrap.empirical.95.genome.comparisons, file.path(
  ANCESTRY.OUTPUT.DIR, ancestry.output.filename(
    "ancestry.bootstrap.empirical.95.genome.comparisons", "rds"
    )
  ))
saveRDS(ancestry.bootstrap.empirical.bonferroni.chromosome.comparisons,
  file.path(
    ANCESTRY.OUTPUT.DIR,
    ancestry.output.filename(
      "ancestry.bootstrap.empirical.bonferroni.chromosome.comparisons", "rds"
      )
    )
  )
saveRDS(ancestry.bootstrap.empirical.bonferroni.genome.comparisons, file.path(
  ANCESTRY.OUTPUT.DIR, ancestry.output.filename(
    "ancestry.bootstrap.empirical.bonferroni.genome.comparisons", "rds"
    )
  ))
saveRDS(ancestry.bootstrap.simulation.95.comparisons, file.path(
  ANCESTRY.OUTPUT.DIR, ancestry.output.filename(
    "ancestry.bootstrap.simulation.95.comparisons", "rds"
    )
  ))
saveRDS(ancestry.bootstrap.simulation.bonferroni.comparisons, file.path(
  ANCESTRY.OUTPUT.DIR, ancestry.output.filename(
    "ancestry.bootstrap.simulation.bonferroni.comparisons", "rds"
    )
  ))
print(ancestry.bootstrap.tcd.1kg.bar)
print(ancestry.bootstrap.all.datatypes.adx.asw.bar)
print(ancestry.bootstrap.tcd.1kg.histograms$chr1)
print(ancestry.bootstrap.empirical.bonferroni.genome.comparisons)
print(ancestry.bootstrap.simulation.bonferroni.comparisons)

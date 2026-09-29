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
library(tidyverse)
library(ggh4x)
library(ggfx)
library(nanoparquet)


SIM.TC.DATA.DIR <- "~/scratch/OOA_NAAdmixture_2T12Consistent/stats"
SIMDOWN.TC.DATA.DIR <-
  "~/scratch/OOA_NAAdmixture_2T12ConsistentOnekgDownsample/stats"
SIM.LG.DATA.DIR <- "~/scratch/OOA_NAAdmixture_largeGrowth/stats"
SIMDOWN.LG.DATA.DIR <-
  "~/scratch/OOA_NAAdmixture_largeGrowthOnekgDownsample/stats"
EMPIRICAL.DATA.DIR <- "~/scratch/OOA_NAAdmixture_1kG/stats"
OUTPUT.DIR <- "/home1/karatas/proj/OOA_NAAdmixture_data"
CHROMOSOMES <- as.character(1:22)
SELECTED.CHROMOSOMES <- c("1", "10", "20")
SIMULATION.K <- 2L
EMPIRICAL.K <- 2L
RANDOM.SEED <- 123L
DOWNSAMPLE.SIZE <- 50L
BOOTSTRAP.REPLICATES <- 1000L
ANCESTRY.BOOTSTRAP.CHROMOSOME.FAMILY.SIZE <- 308L
ANCESTRY.BOOTSTRAP.GENOME.ASW.FAMILY.SIZE <- 176L
HISTOGRAM.BREAKS <- seq(0, 1, by = 0.05)
PLOT.EMPIRICAL.METHOD <- "ADMIXTURE"
PLOT.SAMPLE.SET <- "downsampled"
PLOT.BASE.SIZE <- 24
CATEGORICAL.BAR.DODGE <- 0.9
CATEGORICAL.BAR.WIDTH <- 0.8
CATEGORICAL.BAR.LINEWIDTH <- 1
DENSE.BAR.WIDTH.MULTIPLIER <- 0.8
DENSE.BAR.LINEWIDTH <- 0.75
SOURCE.LEVELS <- c(
  "Simulation_2T12Consistent", "Simulation_2T12Consistent_simDown",
  "Simulation_largeGrowth", "Simulation_largeGrowth_simDown", "Empirical"
  )
SOURCE.LABELS <- setNames(
  c("T.C.", "T.C.D.", "L.G.", "L.G.D.", "Emp."), SOURCE.LEVELS
  )
PLOT.STYLES <- list(
  colors = c(
    Simulation_2T12Consistent = "#9A83CE",
    Simulation_2T12Consistent_simDown = "#6F55B5",
    Simulation_largeGrowth = "#32146F",
    Simulation_largeGrowth_simDown = "#4B1FA8", Empirical = "#B83264"
    ),
  labels = SOURCE.LABELS,
  empirical.colors = c(ADMIXTURE = "#B83264", fastStructure = "#B9584A"),
  contrast.colors = c(
    `TC-Emp` = "#9A83CE", `TCD-Emp` = "#6F55B5",
    `LG-Emp` = "#32146F", `LGD-Emp` = "#4B1FA8",
    `TC-TCD` = "#BDBDBD", `LG-LGD` = "#737373",
    `TC-LG` = "#000000"
    ),
  contrast.labels = c(
    `TC-Emp` = "T.C. - ASW", `TCD-Emp` = "T.C.D. - ASW",
    `LG-Emp` = "L.G. - ASW", `LGD-Emp` = "L.G.D. - ASW",
    `TC-TCD` = "T.C. - T.C.D.", `LG-LGD` = "L.G. - L.G.D.",
    `TC-LG` = "T.C. - L.G."
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
      (data.type %in% SOURCE.LEVELS[1:4] & pop == "ADX") |
        (data.type == "Empirical" & pop == "ASW")
      ) %>%
    mutate(data.type = factor(data.type, levels = SOURCE.LEVELS))
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
      (data.type %in% SOURCE.LEVELS[c(1, 3)] & method == "tspop") |
        (data.type %in% SOURCE.LEVELS[c(2, 4)] &
          method == PLOT.EMPIRICAL.METHOD)
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
      method == PLOT.EMPIRICAL.METHOD
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
  selected <- select.bootstrap.ancestry.ids(data, DOWNSAMPLE.SIZE,
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
      method == PLOT.EMPIRICAL.METHOD
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
      chrom %in% SELECTED.CHROMOSOMES
      ) %>%
    mutate(data.type = factor(
      as.character(data.type),
      levels = order.active.levels(data.type, SOURCE.LEVELS)
      ))
  genome <- data %>%
    filter(data.type == "Empirical", chrom == "all", stat %in% plotted$stat)
  dodge <- position_dodge(width = CATEGORICAL.BAR.DODGE)
  plot <- ggplot(plotted, aes(chrom, mean, fill = data.type)) +
    geom_rect(
      data = genome, aes(ymin = lower, ymax = upper),
      xmin = -Inf, xmax = Inf, inherit.aes = FALSE,
      fill = PLOT.STYLES$empirical.colors[[PLOT.EMPIRICAL.METHOD]],
      alpha = 0.15
      ) +
    geom_col(
      position = dodge, width = CATEGORICAL.BAR.WIDTH,
      linewidth = CATEGORICAL.BAR.LINEWIDTH, color = "black"
      ) +
    geom_errorbar(
      aes(ymin = lower, ymax = upper), position = dodge, width = 0,
      linewidth = CATEGORICAL.BAR.LINEWIDTH, na.rm = TRUE
      ) +
    with_outer_glow(
      geom_hline(
        data = genome, aes(yintercept = mean), linetype = "longdash",
        color = PLOT.STYLES$empirical.colors[[PLOT.EMPIRICAL.METHOD]],
        linewidth = 1
        ),
      colour = "black", sigma = 0, expand = 3
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
    scale_fill_manual(values = PLOT.STYLES$colors,
      labels = PLOT.STYLES$labels) +
    labs(title = title, x = "Chromosome", y = "African ancestry", fill = NULL) +
    theme_bw(base_size = PLOT.BASE.SIZE) +
    theme(legend.position = "top", panel.grid.minor = element_blank())
  return(plot)
  }


# construct chromosome-1 simulation and genome-wide asw histograms.
make.bootstrap.ancestry.histogram.plot <- function(data, data.types, title) {
  plotted <- data %>%
    filter(
      as.character(data.type) %in% data.types,
      (as.character(chrom) == "all" & data.type == "Empirical") |
        (as.character(chrom) == "1" & data.type != "Empirical")
      ) %>%
    mutate(data.type = factor(
      as.character(data.type),
      levels = order.active.levels(data.type, SOURCE.LEVELS)
      ))
  dodge <- position_dodge(width = diff(HISTOGRAM.BREAKS)[1])
  plot <- ggplot(plotted, aes(xmid, mean, fill = data.type)) +
    geom_col(
      position = dodge,
      width = diff(HISTOGRAM.BREAKS)[1] * DENSE.BAR.WIDTH.MULTIPLIER,
      color = "black", linewidth = DENSE.BAR.LINEWIDTH
      ) +
    geom_errorbar(
      aes(ymin = lower, ymax = upper), position = dodge, width = 0,
      linewidth = DENSE.BAR.LINEWIDTH, na.rm = TRUE
      ) +
    scale_fill_manual(values = PLOT.STYLES$colors,
      labels = PLOT.STYLES$labels) +
    labs(
      title = title, x = "African ancestry", y = "Fraction of individuals",
      fill = NULL
      ) +
    theme_bw(base_size = PLOT.BASE.SIZE) +
    theme(legend.position = "top", panel.grid.minor = element_blank())
  return(plot)
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
    reference.scope, family.size, bootstrap.replicates
  ) {
  table <- map_dfr(chromosomes, function(chromosome) {
    map_dfr(seq_len(nrow(contrasts)), function(index) {
      contrast <- contrasts[index, ]
      left.values <- simulation %>%
        filter(source == contrast$left.source, chrom == chromosome) %>%
        arrange(rep) %>% pull(all_of(statistic))
      if (contrast$right.source == "Emp") {
        empirical.values <- empirical %>%
          filter(chrom == reference.chrom(chromosome)) %>% pull(afr.q)
        bootstrap <- bootstrap.ancestry.empirical.difference(
          left.values, empirical.values, statistic, bootstrap.replicates
          )
        } else {
        right.values <- simulation %>%
          filter(source == contrast$right.source, chrom == chromosome) %>%
          arrange(rep) %>% pull(all_of(statistic))
        bootstrap <- bootstrap.ancestry.simulation.difference(
          left.values, right.values, contrast$paired, bootstrap.replicates
          )
        }
      return(bind_cols(
        tibble(
          statistic = statistic, contrast = contrast$contrast,
          chromosome = chromosome, reference.scope = reference.scope
          ),
        summarize.ancestry.bootstrap.difference(
          bootstrap$draws, bootstrap$difference, family.size
          )
        ))
      })
    })
  return(table)
  }


# make deterministic chromosome and genome-wide asw comparison tables.
make.ancestry.bootstrap.comparison.tables <- function(
    summary.data, individual.data, empirical.method = PLOT.EMPIRICAL.METHOD,
    sample.set.input = PLOT.SAMPLE.SET, chromosomes = CHROMOSOMES,
    bootstrap.replicates = BOOTSTRAP.REPLICATES, seed = RANDOM.SEED
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
  set.seed(seed)
  chromosome.tables <- map_dfr(c("mean", "sd"), function(statistic) {
    make.ancestry.bootstrap.comparison.table(
      data$simulation, data$empirical, ANCESTRY.BOOTSTRAP.CONTRASTS,
      chromosomes, statistic, identity, "chromosome",
      ANCESTRY.BOOTSTRAP.CHROMOSOME.FAMILY.SIZE, bootstrap.replicates
      )
    })
  genome.tables <- map_dfr(c("mean", "sd"), function(statistic) {
    make.ancestry.bootstrap.comparison.table(
      data$simulation, data$empirical, genome.contrasts, chromosomes,
      statistic, function(chromosome) "all", "genome",
      ANCESTRY.BOOTSTRAP.GENOME.ASW.FAMILY.SIZE, bootstrap.replicates
      )
    })
  if (nrow(chromosome.tables) != ANCESTRY.BOOTSTRAP.CHROMOSOME.FAMILY.SIZE ||
      nrow(genome.tables) != ANCESTRY.BOOTSTRAP.GENOME.ASW.FAMILY.SIZE) {
    stop("Bootstrap ancestry comparison family size does not match its table")
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
      "ancestry.statistics.bootstrap.chromosome.comparisons.csv"
      )
    )
  readr::write_csv(
    tables$chromosome.vs.genome.asw,
    file.path(
      output.directory,
      "ancestry.statistics.bootstrap.chromosome_vs_genome_asw.csv"
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
      chromosome.levels = CHROMOSOMES
      ),
    simulation = list(
      contrasts = c("TC-TCD", "LG-LGD", "TC-LG"),
      chromosome.levels = CHROMOSOMES
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
      contrast.color = unname(PLOT.STYLES$contrast.colors[as.character(
        contrast
        )]),
      facet.label = factor(NA_character_)
      ) %>%
    group_by(chromosome) %>%
    mutate(
      plot.x = as.numeric(chromosome) +
        (as.numeric(contrast) - (n_distinct(contrast) + 1) / 2) *
          CATEGORICAL.BAR.DODGE / n_distinct(contrast)
      ) %>%
    ungroup() %>%
    arrange(statistic, chromosome, contrast)
  return(plotted)
  }


# construct a faceted bootstrap contrast plot from prepared comparison data.
make.ancestry.bootstrap.contrast.plot <- function(
    data, title, facet.statistics = TRUE
  ) {
  plot <- ggplot(
    data,
    aes(plot.x, difference, color = contrast, group = contrast)
    ) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_errorbar(
      aes(ymin = ci.lower, ymax = ci.upper),
      width = CATEGORICAL.BAR.WIDTH * 0.2,
      linewidth = CATEGORICAL.BAR.LINEWIDTH
      ) +
    geom_point(size = 2.2) +
    geom_point(
      data = filter(data, significant),
      aes(x = plot.x, y = difference, fill = contrast, group = contrast),
      shape = 21, color = "red", size = 2.2,
      stroke = CATEGORICAL.BAR.LINEWIDTH, inherit.aes = FALSE
      ) +
    scale_color_manual(
      values = PLOT.STYLES$contrast.colors,
      labels = PLOT.STYLES$contrast.labels
      ) +
    scale_fill_manual(
      values = PLOT.STYLES$contrast.colors, guide = "none"
      ) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
    labs(
      x = "Chromosome",
      y = "Difference", color = NULL,
      title = title
      ) +
    scale_x_continuous(
      breaks = seq_along(levels(data$chromosome)),
      labels = levels(data$chromosome)
      ) +
    theme_bw(base_size = PLOT.BASE.SIZE) +
    theme(
      legend.position = "top",
      legend.title = element_blank(),
      strip.background = element_rect(fill = "grey92"),
      panel.grid.minor = element_blank()
      )
    
  if (facet.statistics) {
    if (all(is.na(data$facet.label))) {
      plot <- plot + facet_grid(
        rows = vars(statistic), scales = "free_y",
        labeller = as_labeller(c(mean = "Mean", sd = "SD"))
        )
      } else {
      plot <- plot + facet_wrap(
        vars(facet.label), ncol = 1, scales = "free_y",
        strip.position = "right"
        )
      }
    }
  return(plot)
  }


# analysis ----


# read tc/lg truth, configured tcd/lgd inference, and empirical inference.
ancestry.sim.tc.tspop.data <- read.ancestry.family(
  SIM.TC.DATA.DIR, "ancestry.chr{chrom}.parquet", CHROMOSOMES,
  "Simulation_2T12Consistent", "tspop", 0L, "tspop"
  )
ancestry.simDown.tc.inference.data <- read.ancestry.family(
  SIMDOWN.TC.DATA.DIR, ancestry.inference.file.family(PLOT.EMPIRICAL.METHOD),
  CHROMOSOMES, "Simulation_2T12Consistent_simDown", PLOT.EMPIRICAL.METHOD,
  SIMULATION.K, PLOT.EMPIRICAL.METHOD
  )
ancestry.sim.lg.tspop.data <- read.ancestry.family(
  SIM.LG.DATA.DIR, "ancestry.chr{chrom}.parquet", CHROMOSOMES,
  "Simulation_largeGrowth", "tspop", 0L, "tspop"
  )
ancestry.simDown.lg.inference.data <- read.ancestry.family(
  SIMDOWN.LG.DATA.DIR, ancestry.inference.file.family(PLOT.EMPIRICAL.METHOD),
  CHROMOSOMES, "Simulation_largeGrowth_simDown", PLOT.EMPIRICAL.METHOD,
  SIMULATION.K, PLOT.EMPIRICAL.METHOD
  )
ancestry.empirical.file.family <- ancestry.inference.file.family(
  PLOT.EMPIRICAL.METHOD
  )
ancestry.empirical.genome.file <- str_replace(
  ancestry.empirical.file.family, "\\.chr\\{chrom\\}", ""
  )
ancestry.empirical.inference.data <- read.ancestry.family(
  EMPIRICAL.DATA.DIR, ancestry.empirical.file.family, CHROMOSOMES, "Empirical",
  PLOT.EMPIRICAL.METHOD, EMPIRICAL.K, "Empirical", include.genome = TRUE,
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
  ancestry.individual.data, DOWNSAMPLE.SIZE, "sample_id",
  c("data.type", "rep", "chrom"), RANDOM.SEED
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
  empirical.method = PLOT.EMPIRICAL.METHOD,
  sample.set.input = PLOT.SAMPLE.SET, chromosomes = CHROMOSOMES,
  bootstrap.replicates = BOOTSTRAP.REPLICATES, seed = RANDOM.SEED
  )
write.ancestry.bootstrap.comparison.tables(
  ancestry.bootstrap.comparison.tables, OUTPUT.DIR
  )


# construct the four active bottom bootstrap figures.
ancestry.bootstrap.summary <- summarize.bootstrap.ancestry(
  ancestry.individual.data, DOWNSAMPLE.SIZE, RANDOM.SEED,
  BOOTSTRAP.REPLICATES
  )
ancestry.bootstrap.histograms <- summarize.bootstrap.histograms(
  ancestry.individual.data, HISTOGRAM.BREAKS, RANDOM.SEED,
  BOOTSTRAP.REPLICATES
  )
ancestry.bootstrap.tcd.1kg.bar <- make.bootstrap.ancestry.bar.plot(
  ancestry.bootstrap.summary,
  c("Simulation_2T12Consistent_simDown", "Empirical"),
  "African ancestry: TCD and ASW"
  )
ancestry.bootstrap.all.datatypes.adx.asw.bar <-
  make.bootstrap.ancestry.bar.plot(
  ancestry.bootstrap.summary, SOURCE.LEVELS,
  "African ancestry: all ADX sources and ASW"
  )
ancestry.bootstrap.tcd.1kg.histogram <- make.bootstrap.ancestry.histogram.plot(
  ancestry.bootstrap.histograms,
  c("Simulation_2T12Consistent_simDown", "Empirical"),
  "Chromosome 1 simulations and genome-wide ASW: TCD and ASW"
  )
ancestry.bootstrap.all.datatypes.adx.asw.histogram <-
  make.bootstrap.ancestry.histogram.plot(
    ancestry.bootstrap.histograms, SOURCE.LEVELS,
    "Chromosome 1 simulations and genome-wide ASW: all ADX sources and ASW"
    )
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
dir.create(OUTPUT.DIR, recursive = TRUE, showWarnings = FALSE)
saveRDS(ancestry.bootstrap.tcd.1kg.bar, file.path(
  OUTPUT.DIR, "ancestry.bootstrap.tcd.1kg.bar.rds"
  ))
saveRDS(ancestry.bootstrap.all.datatypes.adx.asw.bar, file.path(
  OUTPUT.DIR, "ancestry.bootstrap.all.datatypes.adx.asw.bar.rds"
  ))
saveRDS(ancestry.bootstrap.tcd.1kg.histogram, file.path(
  OUTPUT.DIR, "ancestry.bootstrap.tcd.1kg.histogram.rds"
  ))
saveRDS(ancestry.bootstrap.all.datatypes.adx.asw.histogram, file.path(
  OUTPUT.DIR, "ancestry.bootstrap.all.datatypes.adx.asw.histogram.rds"
  ))
saveRDS(ancestry.bootstrap.empirical.95.chromosome.comparisons, file.path(
  OUTPUT.DIR, "ancestry.bootstrap.empirical.95.chromosome.comparisons.rds"
  ))
saveRDS(ancestry.bootstrap.empirical.95.genome.comparisons, file.path(
  OUTPUT.DIR, "ancestry.bootstrap.empirical.95.genome.comparisons.rds"
  ))
saveRDS(ancestry.bootstrap.empirical.bonferroni.chromosome.comparisons,
  file.path(
    OUTPUT.DIR,
    "ancestry.bootstrap.empirical.bonferroni.chromosome.comparisons.rds"
    )
  )
saveRDS(ancestry.bootstrap.empirical.bonferroni.genome.comparisons, file.path(
  OUTPUT.DIR, "ancestry.bootstrap.empirical.bonferroni.genome.comparisons.rds"
  ))
saveRDS(ancestry.bootstrap.simulation.95.comparisons, file.path(
  OUTPUT.DIR, "ancestry.bootstrap.simulation.95.comparisons.rds"
  ))
saveRDS(ancestry.bootstrap.simulation.bonferroni.comparisons, file.path(
  OUTPUT.DIR, "ancestry.bootstrap.simulation.bonferroni.comparisons.rds"
  ))
print(ancestry.bootstrap.tcd.1kg.bar)
print(ancestry.bootstrap.all.datatypes.adx.asw.bar)
print(ancestry.bootstrap.tcd.1kg.histogram)
print(ancestry.bootstrap.all.datatypes.adx.asw.histogram)
print(ancestry.bootstrap.empirical.95.chromosome.comparisons)
print(ancestry.bootstrap.empirical.95.genome.comparisons)
print(ancestry.bootstrap.empirical.bonferroni.chromosome.comparisons)
print(ancestry.bootstrap.empirical.bonferroni.genome.comparisons)
print(ancestry.bootstrap.simulation.95.comparisons)
print(ancestry.bootstrap.simulation.bonferroni.comparisons)

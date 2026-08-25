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

SIM.SMALL.DATA.DIR <- "~/scratch/OOA_NAAdmixture_small/stats"
SIM.LARGE.DATA.DIR <- "~/scratch/OOA_NAAdmixture_large/stats"
EMPIRICAL.DATA.DIR <- "~/scratch/OOA_NAAdmixture_1kG_cp/stats"
CHROMOSOME.LENGTHS.PATH <- "~/proj/1000GenomeNYGC_hg38_karatas/ONEKG_chr_lens.tsv"
CHROMOSOMES <- c("1", "5", "10", "14", "18", "22")
SIMULATION.K <- 2
EMPIRICAL.K <- 2
RANDOM.SEED <- 123
DOWNSAMPLE.SIZE <- 50
BOOTSTRAP.REPLICATES <- 1000
HISTOGRAM.BREAKS <- seq(0, 1, by = 0.05)
ADMIXED.ROLES <- c("ADX", "ASW")
ANCESTRY.COMPONENT.COLORS <- c(
  component_1_q = "#0072B2", component_2_q = "#D55E00"
)
PLOT.STYLES <- list(
  colors = c(
    Simulation_small.full = "#00AEDB",
    Simulation_small.downsampled = "#75D5F0",
    Simulation_large.full = "#007E9F",
    Simulation_large.downsampled = "#6BAFC1",
    Empirical.full = "#B83264"
  ),
  labels = c(
    Simulation_small.full = "Simulation small",
    Simulation_small.downsampled = "Downsampled small",
    Simulation_large.full = "Simulation large",
    Simulation_large.downsampled = "Downsampled large",
    Empirical.full = "Empirical"
  ),
  shapes = c(full = 21, downsampled = 24),
  linetypes = c(full = "solid", downsampled = "dashed"),
  simulation.labels = c(
    tspop = "tspop", ADMIXTURE = "ADMIXTURE",
    fastStructure = "fastStructure"
  ),
  empirical.colors = c(
    ADMIXTURE = "#B83264", fastStructure = "#B9584A"
  )
)


# internal functions ----


# standardize one ancestry table to the shared analysis schema
normalize.ancestry.table <- function(
    data, data.type.input, method, simulation.source.input, k,
    population.mapping = NULL, sample.id.column = "sample_id"
) {
  # add shared identifiers, metadata, and optional population annotations
  data$chrom <- as.character(data$chrom)
  data$sample_id <- as.character(data[[sample.id.column]])
  data$rep <- if ("rep" %in% names(data)) as.numeric(data$rep) else 0
  data$data.type <- data.type.input
  data$method <- method
  data$simulation.source <- simulation.source.input
  if (!is.null(population.mapping)) {
    data <- left_join(data, population.mapping, by = "sample_id")
  }
  if (!"pop" %in% names(data)) data$pop <- NA_character_
  if (!"role" %in% names(data)) data$role <- data$pop
  if (!"vcf_sample_id" %in% names(data)) {
    data$vcf_sample_id <- data$sample_id
  }
  if (!"span" %in% names(data)) data$span <- NA_real_

  # normalize tspop truth or retain only the requested inference K
  if (method == "tspop") {
    data$component_1_q <- as.numeric(data$afr_tspop)
    data$component_2_q <- as.numeric(data$eur_tspop)
    data <- select(data, -afr_tspop, -eur_tspop)
    data$k <- 0
  } else {
    components <- paste0("component_", seq_len(k), "_q")
    if (!"k" %in% names(data)) data$k <- k
    data <- filter(data, k == !!k)
  }
  # pad unused components and enforce numeric ancestry columns
  for (component in paste0("component_", 3:5, "_q")) {
    if (!component %in% names(data)) data[[component]] <- NA_real_
  }
  data <- data %>% mutate(
    component_1_q = as.numeric(component_1_q),
    component_2_q = as.numeric(component_2_q), span = as.numeric(span)
  )

  return(data)
}


# read chromosome or genome files for one ancestry source and normalize them
read.ancestry.family <- function(
    data.directory, file.family, chromosomes, data.type.input, method, k,
    simulation.source.input, population.mapping = NULL,
    sample.id.column = "sample_id", include.genome = FALSE,
    genome.file.family = NULL
) {
  # expand the filename pattern into concrete chromosome paths
  paths <- file.path(
    path.expand(data.directory),
    vapply(chromosomes, function(x) {
      return(gsub("\\{chrom\\}", x, file.family))
    }, character(1))
  )
  # append the whole-genome file when requested for empirical data
  if (include.genome) {
    chromosomes <- c(chromosomes, "all")
    paths <- c(paths, file.path(
      path.expand(data.directory), genome.file.family
    ))
  }
  # read, label, normalize, and combine every requested table
  data <- map2_dfr(paths, chromosomes, function(path, chrom) {
    table <- read_parquet(path)
    table$chrom <- chrom
    return(normalize.ancestry.table(
      table, data.type.input, method, simulation.source.input, k,
      population.mapping, sample.id.column
    ))
  })

  return(data)
}


# orient K=2 components so afr.q consistently represents African ancestry
orient.ancestry.components <- function(
    data, admixed.roles, grouping.columns
) {
  # retain tspop truth and the configured K for inference methods
  filtered <- filter(
    data, method == "tspop" | k == SIMULATION.K | k == EMPIRICAL.K
  )
  # determine the African component within each analysis group
  component.map <- filtered %>%
    filter(role %in% admixed.roles) %>%
    group_by(across(all_of(grouping.columns))) %>%
    summarise(
      component.1 = mean(component_1_q, na.rm = TRUE),
      component.2 = mean(component_2_q, na.rm = TRUE), .groups = "drop"
    ) %>%
    mutate(afr.component = if_else(
      component.1 >= component.2, "component_1_q", "component_2_q"
    ))
  # apply the orientation and discard the temporary component map
  data <- filtered %>%
    left_join(component.map, by = grouping.columns) %>%
    mutate(
      afr.q = if_else(
        afr.component == "component_1_q", component_1_q, component_2_q
      ),
      eur.q = if_else(
        afr.component == "component_1_q", component_2_q, component_1_q
      )
    ) %>%
    select(-component.1, -component.2, -afr.component)

  return(data)
}


# select reproducible IDs per simulation size, replicate, and chromosome
select.downsample.ids <- function(
    data, downsample.size, sample.id.column, grouping.columns, seed
) {
  # use tspop identifiers as the candidates shared across all methods
  candidates <- data %>%
    filter(data.type != "Empirical", method == "tspop") %>%
    distinct(across(all_of(c(grouping.columns, sample.id.column))))
  # Verify every simulation group can supply the fixed sample size.
  sizes <- candidates %>%
    count(across(all_of(grouping.columns)), name = "available")
  if (any(sizes$available < downsample.size)) {
    stop(paste0(
      "A simulation group contains fewer than ", downsample.size,
      " candidates"
    ))
  }
  # draw reproducibly while allowing selections to vary by replicate
  set.seed(seed)
  selected.ids <- candidates %>%
    group_by(across(all_of(grouping.columns))) %>%
    mutate(.draw.group = cur_group_id(), .random.order = runif(n())) %>%
    arrange(.random.order, .by_group = TRUE) %>%
    group_modify(function(group, key) {
      start <- ((first(group$.draw.group) - 1) %% nrow(group)) + 1
      indices <- ((start - 1 + seq_len(downsample.size) - 1) %%
        nrow(group)) + 1
      return(slice(group, indices))
    }) %>%
    ungroup() %>%
    select(-.draw.group, -.random.order) %>%
    select(all_of(c(grouping.columns, sample.id.column))) %>%
    arrange(across(all_of(grouping.columns)), .data[[sample.id.column]])

  return(selected.ids)
}


# add full and downsampled rows using the shared selected identifiers
apply.downsample.ids <- function(
    data, selected.ids, sample.id.column, grouping.columns
) {
  # reject duplicate selections that would multiply joined observations
  duplicates <- selected.ids %>%
    count(across(all_of(c(grouping.columns, sample.id.column)))) %>%
    filter(n != 1)
  if (nrow(duplicates)) stop("Selected sample IDs are not unique")
  # retain all rows as full and duplicate selected simulation rows only
  full <- mutate(data, sample.set = "full")
  downsampled <- data %>%
    filter(data.type != "Empirical") %>%
    inner_join(selected.ids, by = c(grouping.columns, sample.id.column)) %>%
    mutate(sample.set = "downsampled")
  data <- bind_rows(full, downsampled)

  return(data)
}


# calculate a percentile bootstrap interval for a supplied statistic.
bootstrap.interval <- function(values, statistic, replicates, seed) {
  # return undefined bounds when the sample or bootstrap is too small.
  if (length(values) < 2 || replicates < 2) return(c(NA_real_, NA_real_))
  # resample deterministically and extract the central 95% interval.
  set.seed(seed)
  estimates <- replicate(replicates, {
    statistic(sample(values, length(values), replace = TRUE))
  })
  interval <- as.numeric(
    quantile(estimates, c(0.025, 0.975), names = FALSE)
  )

  return(interval)
}


# summarize ancestry and add empirical bootstrap uncertainty intervals
summarize.ancestry <- function(
    data, grouping.columns, bootstrap.replicates, bootstrap.seed,
    admixed.roles = "ADX"
) {
  # compute descriptive statistics within each requested analysis group
  result <- data %>%
    filter(role %in% admixed.roles, !is.na(afr.q)) %>%
    group_by(across(all_of(grouping.columns))) %>%
    group_modify(function(group, key) {
      x <- group$afr.q
      empirical <- key$data.type[[1]] == "Empirical"
      mean.ci <- if (empirical) bootstrap.interval(
        x, mean, bootstrap.replicates, bootstrap.seed
      ) else c(NA_real_, NA_real_)
      sd.ci <- if (empirical) bootstrap.interval(
        x, sd, bootstrap.replicates, bootstrap.seed + 1
      ) else c(NA_real_, NA_real_)
      return(tibble(
        mean = mean(x), sd = sd(x), median = median(x),
        q25 = quantile(x, 0.25, names = FALSE),
        q75 = quantile(x, 0.75, names = FALSE), n = length(x),
        mean.boot = if (empirical) mean(x) else NA_real_,
        mean.boot.lower = mean.ci[1], mean.boot.upper = mean.ci[2],
        sd.boot = if (empirical) sd(x) else NA_real_,
        sd.boot.lower = sd.ci[1], sd.boot.upper = sd.ci[2]
      ))
    }) %>%
    ungroup()
  # reshape downsampled statistics for side-by-side inspection
  downsampled <- result %>%
    filter(sample.set == "downsampled") %>%
    select(-sample.set, -contains("boot")) %>%
    rename(
      mean.rand.downsample = mean, sd.rand.downsample = sd,
      median.rand.downsample = median, q25.rand.downsample = q25,
      q75.rand.downsample = q75, n.rand.downsample = n
    )
  # join downsampled values onto their corresponding summary groups
  result <- left_join(result, downsampled,
    by = setdiff(grouping.columns, "sample.set")
  )

  return(result)
}


# convert ancestry values into replicate-aggregated histogram bins
summarize.histograms <- function(
    data, breaks, chromosomes, admixed.roles = "ADX"
) {
  # keep valid admixed values and empirical-only whole-genome rows
  histogram.data <- data %>%
    filter(
      role %in% admixed.roles, !is.na(afr.q), between(afr.q, 0, 1),
      chrom %in% c(chromosomes, "all"),
      !(chrom == "all" & data.type != "Empirical")
    ) %>%
    # calculate normalized bin counts independently for each replicate
    mutate(chrom = factor(chrom, levels = c(chromosomes, "all"))) %>%
    group_by(chrom, rep, data.type, simulation.source, method, sample.set) %>%
    group_modify(function(group, key) {
      h <- hist(group$afr.q, breaks = breaks, plot = FALSE)
      return(tibble(
        xmin = head(h$breaks, -1), xmax = tail(h$breaks, -1),
        xmid = h$mids, fraction = h$counts / sum(h$counts)
      ))
    }) %>%
    ungroup() %>%
    # aggregate replicate fractions and calculate simulation error bounds
    group_by(
      chrom, data.type, simulation.source, method, sample.set,
      xmin, xmax, xmid
    ) %>%
    summarise(
      mean.frac = mean(fraction),
      sd.frac = if (n() > 1) sd(fraction) else NA_real_,
      n.rep = n_distinct(rep), .groups = "drop"
    ) %>%
    mutate(
      ymin = if_else(data.type == "Empirical", NA_real_,
        pmax(0, mean.frac - 2 * sd.frac)
      ),
      ymax = if_else(data.type == "Empirical", NA_real_,
        mean.frac + 2 * sd.frac
      )
    )

  return(histogram.data)
}


# define the four allowed simulation and empirical method combinations
plot.variant.specifications <- function() {
  # exclude cross-method inference combinations by construction
  specifications <- tribble(
    ~key, ~simulation.source, ~empirical.method,
    "tspop.admixture", "tspop", "ADMIXTURE",
    "admixture.admixture", "ADMIXTURE", "ADMIXTURE",
    "tspop.fastStructure", "tspop", "fastStructure",
    "fastStructure.fastStructure", "fastStructure", "fastStructure"
  )

  return(specifications)
}


# select the five data series used by one plot combination
prepare.plot.data <- function(
    summary.data, simulation.source.input, empirical.method, chromosomes
) {
  # pair one simulation source with one complete empirical method
  plot.data <- summary.data %>%
    filter(
      (data.type != "Empirical" &
        .data$simulation.source == simulation.source.input &
        chrom %in% chromosomes) |
        (data.type == "Empirical" & method == empirical.method &
          sample.set == "full" & chrom %in% c(chromosomes, "all"))
    ) %>%
    mutate(
      chrom = factor(chrom, levels = c(chromosomes, "all")),
      series = paste(data.type, sample.set, sep = ".")
    )

  return(plot.data)
}


# build a chromosome plot for either mean or standard deviation
make.stat.by.chrom.plot <- function(
    summary.data, simulation.source.input, empirical.method, chromosomes,
    styles, statistic, y.label
) {
  # separate simulation distributions from empirical reference values
  data <- prepare.plot.data(
    summary.data, simulation.source.input, empirical.method, chromosomes
  )
  simulation <- filter(data, data.type != "Empirical")
  empirical <- filter(data, data.type == "Empirical")
  genome <- filter(empirical, chrom == "all")
  empirical.chrom <- filter(empirical, chrom != "all")
  lower <- paste0(statistic, ".boot.lower")
  upper <- paste0(statistic, ".boot.upper")
  estimate <- paste0(statistic, ".boot")
  color <- styles$empirical.colors[[empirical.method]]
  # draw simulation boxes with empirical chromosome and genome uncertainty
  plot <- ggplot(simulation, aes(chrom, .data[[statistic]], fill = series)) +
    geom_rect(data = genome,
      aes(xmin = -Inf, xmax = Inf, ymin = .data[[lower]],
        ymax = .data[[upper]]),
      inherit.aes = FALSE, fill = color, alpha = 0.1
    ) +
    geom_hline(data = genome, aes(yintercept = .data[[estimate]]),
      color = color, linetype = "dashed"
    ) +
    geom_boxplot(aes(group = interaction(chrom, series)), outlier.shape = NA) +
    geom_errorbar(data = empirical.chrom,
      aes(x = chrom, ymin = .data[[lower]], ymax = .data[[upper]]),
      inherit.aes = FALSE, color = color, width = 0.15
    ) +
    geom_point(data = empirical.chrom, aes(chrom, .data[[estimate]]),
      inherit.aes = FALSE, shape = 23, size = 3, fill = color
    ) +
    scale_x_discrete(limits = chromosomes, drop = FALSE) +
    scale_fill_manual(values = styles$colors, labels = styles$labels) +
    labs(x = "Chromosome", y = y.label, fill = NULL) +
    theme_bw(base_size = 18) +
    theme(legend.position = "top", panel.grid.minor = element_blank())

  return(plot)
}


# build the chromosome-level mean ancestry plot
make.mean.by.chrom.plot <- function(
    summary.data, simulation.source.input, empirical.method, chromosomes, styles
) {
  # delegate construction using the mean statistic and axis label
  plot <- make.stat.by.chrom.plot(
    summary.data, simulation.source.input, empirical.method, chromosomes,
    styles, "mean", "Mean African ancestry"
  )

  return(plot)
}


# build the chromosome-level ancestry standard-deviation plot
make.sd.by.chrom.plot <- function(
    summary.data, simulation.source.input, empirical.method, chromosomes, styles
) {
  # delegate construction using the SD statistic and axis label
  plot <- make.stat.by.chrom.plot(
    summary.data, simulation.source.input, empirical.method, chromosomes,
    styles, "sd", "Standard deviation of African ancestry"
  )

  return(plot)
}


# build vertically faceted mean and standard-deviation chromosome plots
make.mean.sd.plot <- function(
    summary.data, simulation.source.input, empirical.method, chromosomes, styles
) {
  # reshape simulation mean and SD summaries into one plotting table
  data <- prepare.plot.data(
    summary.data, simulation.source.input, empirical.method, chromosomes
  )
  simulation <- data %>%
    filter(data.type != "Empirical") %>%
    pivot_longer(c(mean, sd), names_to = "stat", values_to = "estimate")
  # reshape empirical estimates and bounds to the same statistic key
  empirical <- data %>%
    filter(data.type == "Empirical") %>%
    select(
      chrom, data.type, simulation.source, method, sample.set, series,
      starts_with("mean.boot"), starts_with("sd.boot")
    )
  empirical <- bind_rows(
    transmute(
      empirical, chrom, data.type, simulation.source, method, sample.set,
      series, stat = "mean", estimate = mean.boot,
      lower = mean.boot.lower, upper = mean.boot.upper
    ),
    transmute(
      empirical, chrom, data.type, simulation.source, method, sample.set,
      series, stat = "sd", estimate = sd.boot,
      lower = sd.boot.lower, upper = sd.boot.upper
    )
  )
  color <- styles$empirical.colors[[empirical.method]]
  # draw both statistics with chromosome and genome empirical references
  plot <- ggplot(simulation, aes(chrom, estimate, fill = series)) +
    geom_rect(data = filter(empirical, chrom == "all"),
      aes(xmin = -Inf, xmax = Inf, ymin = lower, ymax = upper),
      inherit.aes = FALSE, fill = color, alpha = 0.1
    ) +
    geom_hline(data = filter(empirical, chrom == "all"),
      aes(yintercept = estimate),
      color = color, linetype = "dashed"
    ) +
    geom_boxplot(aes(group = interaction(chrom, series)), outlier.shape = NA) +
    geom_errorbar(data = filter(empirical, chrom != "all"),
      aes(chrom, ymin = lower, ymax = upper), inherit.aes = FALSE,
      color = color, width = 0.15
    ) +
    geom_point(data = filter(empirical, chrom != "all"),
      aes(chrom, estimate), inherit.aes = FALSE, shape = 23,
      fill = color, size = 3
    ) +
    facet_grid(rows = vars(stat), scales = "free_y") +
    scale_x_discrete(limits = chromosomes, drop = FALSE) +
    scale_fill_manual(values = styles$colors, labels = styles$labels) +
    labs(x = "Chromosome", y = NULL, fill = NULL) +
    theme_bw(base_size = 18)

  return(plot)
}


# relate chromosome length to either mean or SD ancestry summaries
make.length.stat.plot <- function(
    summary.data, chromosome.lengths, simulation.source.input, empirical.method,
    chromosomes, styles, statistic, y.label
) {
  # aggregate replicates and attach the appropriate chromosome length
  data <- prepare.plot.data(
    summary.data, simulation.source.input, empirical.method, chromosomes
  ) %>%
    filter(chrom != "all") %>%
    mutate(chrom = as.character(chrom)) %>%
    group_by(
      chrom, data.type, simulation.source, method, sample.set, series
    ) %>%
    summarise(estimate = median(.data[[statistic]]), .groups = "drop") %>%
    left_join(chromosome.lengths, by = "chrom") %>%
    mutate(chr.len.mb = if_else(
      data.type == "Empirical", chr_len_after_qc, chr_len
    ) / 1e6)
  if (any(is.na(data$chr.len.mb))) {
    stop("Chromosome lengths are unavailable for requested data")
  }
  # draw per-series linear trends and chromosome-level estimates
  plot <- ggplot(data, aes(chr.len.mb, estimate, color = series,
    linetype = sample.set, group = series)) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE) +
    geom_point(aes(shape = sample.set, fill = series), size = 3) +
    scale_color_manual(values = styles$colors, labels = styles$labels) +
    scale_fill_manual(values = styles$colors, labels = styles$labels) +
    scale_shape_manual(values = styles$shapes) +
    scale_linetype_manual(values = styles$linetypes) +
    labs(x = "Chromosome length (Mb)", y = y.label) +
    theme_bw(base_size = 18)

  return(plot)
}


# duild the chromosome-length versus mean-ancestry plot
make.length.mean.plot <- function(
    summary.data, chromosome.lengths, simulation.source.input, empirical.method,
    chromosomes, styles
) {
  # delegate shared length plotting using the mean statistic
  plot <- make.length.stat.plot(
    summary.data, chromosome.lengths, simulation.source.input, empirical.method,
    chromosomes, styles, "mean", "Mean African ancestry"
  )

  return(plot)
}


# build the chromosome-length versus ancestry-SD plot
make.length.sd.plot <- function(
    summary.data, chromosome.lengths, simulation.source.input, empirical.method,
    chromosomes, styles
) {
  # delegate shared length plotting using the standard deviation
  plot <- make.length.stat.plot(
    summary.data, chromosome.lengths, simulation.source.input, empirical.method,
    chromosomes, styles, "sd", "Standard deviation of African ancestry"
  )

  return(plot)
}


# build chromosome histograms with replicate uncertainty for simulations
make.histogram.plot <- function(
    histogram.data, simulation.source.input, empirical.method, chromosomes,
    breaks, styles
) {
  # select one simulation source and its matched empirical method
  data <- histogram.data %>%
    filter(
      (data.type != "Empirical" &
        .data$simulation.source == simulation.source.input &
        chrom %in% chromosomes) |
        (data.type == "Empirical" & method == empirical.method &
          sample.set == "full" & chrom %in% c(chromosomes, "all"))
    ) %>%
    mutate(series = paste(data.type, sample.set, sep = "."))
  error.data <- filter(data, data.type != "Empirical", !is.na(sd.frac))
  dodge <- position_dodge(width = diff(breaks)[1] * 0.95)
  # draw aligned bins, simulation errors, and the empirical all facet
  plot <- ggplot(data, aes(xmid, mean.frac, fill = series, group = series)) +
    geom_col(position = dodge, width = diff(breaks)[1] * 0.95,
      color = "black", linewidth = 0.3) +
    geom_errorbar(data = error.data, aes(ymin = ymin, ymax = ymax),
      position = dodge, width = 0.01) +
    facet_wrap(~ chrom, ncol = 3, drop = TRUE) +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    scale_fill_manual(values = styles$colors, labels = styles$labels) +
    labs(x = "African ancestry",
      y = "Mean fraction of individuals per bin", fill = NULL) +
    theme_bw(base_size = 18)

  return(plot)
}


# build faceted admixture plot for diagnostic inspection
make.diagnostic.admixture.plot <- function(
    individual.data, chromosomes, component.columns, sample.id.column,
    facet.columns, component.colors
) {
  # validate identifiers, component columns, and requested facets
  required <- c("chrom", component.columns, sample.id.column, facet.columns)
  missing <- setdiff(required, names(individual.data))
  if (length(missing)) {
    stop("Diagnostic data are missing: ", paste(missing, collapse = ", "))
  }
  # reshape ancestry components into stacked-bar observations
  data <- individual.data %>%
    filter(chrom %in% chromosomes) %>%
    pivot_longer(all_of(component.columns),
      names_to = "component", values_to = "q")
  # draw one free-width panel for every requested diagnostic group
  plot <- ggplot(data, aes(.data[[sample.id.column]], q, fill = component)) +
    geom_col() +
    facet_wrap(vars(!!!rlang::syms(c(facet.columns, "chrom"))),
      scales = "free_x") +
    scale_fill_manual(values = component.colors) +
    labs(x = NULL, y = "Ancestry proportion", fill = NULL) +
    theme_bw(base_size = 18) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

  return(plot)
}


# read and combine chromosome-labelled fastStructure choose-K diagnostics
read.choose.k.diagnostics <- function(
    data.directory, file.family, chromosomes
) {
  # expand paths and fail before reading an incomplete collection
  paths <- file.path(path.expand(data.directory), vapply(
    chromosomes,
    function(x) return(gsub("\\{chrom\\}", x, file.family)),
    character(1)
  ))
  if (any(!file.exists(paths))) stop("Choose-K files are unavailable")
  # read each table and supply its chromosome when absent
  diagnostics <- map2_dfr(paths, chromosomes, function(path, chrom) {
    data <- read_parquet(path)
    if (!"chrom" %in% names(data)) data$chrom <- chrom
    return(data)
  })

  return(diagnostics)
}


# construct all six plot families from the four allowed combinations
build.plot.families <- function(
    summary.data, histogram.data, chromosome.lengths, chromosomes,
    styles, histogram.breaks = HISTOGRAM.BREAKS
) {
  # pair each family builder with the shared combination specification
  specs <- plot.variant.specifications()
  builders <- list(
    mean.by.chromosome = function(source, method) {
      return(make.mean.by.chrom.plot(
        summary.data, source, method, chromosomes, styles
      ))
    },
    sd.by.chromosome = function(source, method) {
      return(make.sd.by.chrom.plot(
        summary.data, source, method, chromosomes, styles
      ))
    },
    combined.mean.sd = function(source, method) {
      return(make.mean.sd.plot(
        summary.data, source, method, chromosomes, styles
      ))
    },
    length.versus.mean = function(source, method) {
      return(make.length.mean.plot(
        summary.data, chromosome.lengths, source, method,
        chromosomes, styles
      ))
    },
    length.versus.sd = function(source, method) {
      return(make.length.sd.plot(
        summary.data, chromosome.lengths, source, method,
        chromosomes, styles
      ))
    },
    histograms = function(source, method) {
      return(make.histogram.plot(
        histogram.data, source, method, chromosomes,
        histogram.breaks, styles
      ))
    }
  )
  # build and name four ggplot objects within every family
  families <- map(builders, function(builder) {
    plots <- set_names(
      map2(specs$simulation.source, specs$empirical.method, builder),
      specs$key)
    return(plots)
  })

  return(families)
}


# analysis ----


# read simulation ancestry sources
sim.small.tspop.data <- read.ancestry.family(
  SIM.SMALL.DATA.DIR, "ancestry.chr{chrom}.parquet", CHROMOSOMES,
  "Simulation_small", "tspop", 0, "tspop"
)
sim.small.admixture.data <- read.ancestry.family(
  SIM.SMALL.DATA.DIR,
  "ancestry_ADMIXTURE_multik.chr{chrom}.parquet", CHROMOSOMES,
  "Simulation_small", "ADMIXTURE", SIMULATION.K, "ADMIXTURE"
)
sim.small.fastStructure.data <- read.ancestry.family(
  SIM.SMALL.DATA.DIR,
  "ancestry_fastStructure_multik.chr{chrom}.parquet", CHROMOSOMES,
  "Simulation_small", "fastStructure", SIMULATION.K, "fastStructure"
)
sim.large.tspop.data <- read.ancestry.family(
  SIM.LARGE.DATA.DIR, "ancestry.chr{chrom}.parquet", CHROMOSOMES,
  "Simulation_large", "tspop", 0, "tspop"
)
sim.large.admixture.data <- read.ancestry.family(
  SIM.LARGE.DATA.DIR,
  "ancestry_ADMIXTURE_multik.chr{chrom}.parquet", CHROMOSOMES,
  "Simulation_large", "ADMIXTURE", SIMULATION.K, "ADMIXTURE"
)
sim.large.fastStructure.data <- read.ancestry.family(
  SIM.LARGE.DATA.DIR,
  "ancestry_fastStructure_multik.chr{chrom}.parquet", CHROMOSOMES,
  "Simulation_large", "fastStructure", SIMULATION.K, "fastStructure"
)

# read chromosome-level and whole-genome empirical inference
emp.admixture.chromosome.data <- read.ancestry.family(
  EMPIRICAL.DATA.DIR,
  "ancestry_ADMIXTURE_multik.chr{chrom}.parquet", CHROMOSOMES,
  "Empirical", "ADMIXTURE", EMPIRICAL.K, "Empirical"
)
emp.fastStructure.chromosome.data <- read.ancestry.family(
  EMPIRICAL.DATA.DIR,
  "ancestry_fastStructure_multik.chr{chrom}.parquet", CHROMOSOMES,
  "Empirical", "fastStructure", EMPIRICAL.K, "Empirical"
)
emp.admixture.genome.data <- read.ancestry.family(
  EMPIRICAL.DATA.DIR, "ancestry_ADMIXTURE_multik.parquet", "all",
  "Empirical", "ADMIXTURE", EMPIRICAL.K, "Empirical"
)
emp.fastStructure.genome.data <- read.ancestry.family(
  EMPIRICAL.DATA.DIR, "ancestry_fastStructure_multik.parquet", "all",
  "Empirical", "fastStructure", EMPIRICAL.K, "Empirical"
)

# combine, orient, and apply fixed-size downsampling
ancestry.individual.data <- bind_rows(
  sim.small.tspop.data, sim.small.admixture.data,
  sim.small.fastStructure.data, sim.large.tspop.data,
  sim.large.admixture.data, sim.large.fastStructure.data,
  emp.admixture.chromosome.data, emp.fastStructure.chromosome.data,
  emp.admixture.genome.data, emp.fastStructure.genome.data
) %>%
  orient.ancestry.components(
    ADMIXED.ROLES, c("rep", "chrom", "data.type", "method")
  )
downsample.ids <- select.downsample.ids(
  ancestry.individual.data, DOWNSAMPLE.SIZE, "sample_id",
  c("data.type", "rep", "chrom"), RANDOM.SEED
)
ancestry.individual.data <- apply.downsample.ids(
  ancestry.individual.data, downsample.ids, "sample_id",
  c("data.type", "rep", "chrom")
)

# calculate ancestry and histogram summaries
ancestry.summary.data <- summarize.ancestry(
  ancestry.individual.data,
  c(
    "rep", "chrom", "pop", "data.type", "simulation.source",
    "method", "sample.set"
  ),
  BOOTSTRAP.REPLICATES, RANDOM.SEED, ADMIXED.ROLES
)
ancestry.histogram.data <- summarize.histograms(
  ancestry.individual.data, HISTOGRAM.BREAKS, CHROMOSOMES,
  ADMIXED.ROLES
)

# read chromosome lengths and choose-K diagnostics
chromosome.lengths <- readr::read_tsv(
  path.expand(CHROMOSOME.LENGTHS.PATH), show_col_types = FALSE
) %>%
  mutate(chr = str_remove(as.character(chr), "^chr")) %>%
  rename(chrom = chr)
emp.fastStructure.choose.k.chromosome <- read.choose.k.diagnostics(
  EMPIRICAL.DATA.DIR, "fastStructure_chooseK.chr{chrom}.parquet",
  CHROMOSOMES
)
emp.fastStructure.choose.k.genome <- read.choose.k.diagnostics(
  EMPIRICAL.DATA.DIR, "fastStructure_chooseK.parquet", "all"
)
choose.k.frequency.tables <- list(
  chromosome.max.marginal = table(
    emp.fastStructure.choose.k.chromosome$max_marginal_likelihood_k
  ),
  chromosome.model.components = table(
    emp.fastStructure.choose.k.chromosome$model_components_k
  ),
  genome.max.marginal = table(
    emp.fastStructure.choose.k.genome$max_marginal_likelihood_k
  ),
  genome.model.components = table(
    emp.fastStructure.choose.k.genome$model_components_k
  )
)
walk(choose.k.frequency.tables, print)

# create and print diagnostic barplots
simulation.diagnostic.plots <- ancestry.individual.data %>%
  filter(data.type != "Empirical", sample.set == "full") %>%
  split(list(.$data.type, .$method), drop = TRUE) %>%
  map(
    make.diagnostic.admixture.plot,
    chromosomes = CHROMOSOMES,
    component.columns = c("component_1_q", "component_2_q"),
    sample.id.column = "sample_id", facet.columns = c("rep"),
    component.colors = ANCESTRY.COMPONENT.COLORS
  )
empirical.diagnostic.plots <- ancestry.individual.data %>%
  filter(data.type == "Empirical", sample.set == "full") %>%
  split(.$method, drop = TRUE) %>%
  map(
    make.diagnostic.admixture.plot,
    chromosomes = c(CHROMOSOMES, "all"),
    component.columns = c("component_1_q", "component_2_q"),
    sample.id.column = "sample_id", facet.columns = character(),
    component.colors = ANCESTRY.COMPONENT.COLORS
  )
walk(simulation.diagnostic.plots, print)
walk(empirical.diagnostic.plots, print)

# construct the six primary plot-family collections
plot.families <- build.plot.families(
  ancestry.summary.data, ancestry.histogram.data, chromosome.lengths,
  CHROMOSOMES, PLOT.STYLES, HISTOGRAM.BREAKS
)
mean.by.chromosome.plots <- plot.families$mean.by.chromosome
sd.by.chromosome.plots <- plot.families$sd.by.chromosome
combined.mean.sd.plots <- plot.families$combined.mean.sd
length.versus.mean.plots <- plot.families$length.versus.mean
length.versus.sd.plots <- plot.families$length.versus.sd
histogram.plots <- plot.families$histograms

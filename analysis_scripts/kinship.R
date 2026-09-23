# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology
# Mooney Lab
# ___
# kinship.R
# ______________________________________________________________________________

# pattern: Mixed (unavoidable)
# Reason: This analysis script combines pure plot preparation with file I/O.

# set up ----
library(tidyverse)
library(glue)
library(nanoparquet)


SIM.TC.DATA.DIR <- "~/scratch/OOA_NAAdmixture_2T12Consistent/stats"
SIMDOWN.TC.DATA.DIR <- "~/scratch/OOA_NAAdmixture_2T12ConsistentOnekgDownsample/stats"
SIM.LG.DATA.DIR <- "~/scratch/OOA_NAAdmixture_largeGrowth/stats"
SIMDOWN.LG.DATA.DIR <- "~/scratch/OOA_NAAdmixture_largeGrowthOnekgDownsample/stats"
EMPIRICAL.DATA.DIR <- "~/scratch/OOA_NAAdmixture_1kG/stats"
OUTPUT.DIR <- "/home1/karatas/proj/OOA_NAAdmixture_data"
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
  TC.1kG = SOURCE.LEVELS[c(1, 5)],
  TC.TCD = SOURCE.LEVELS[c(1, 2)],
  TCD.1kG = SOURCE.LEVELS[c(2, 5)],
  onlyADX = SOURCE.LEVELS[1:4],
  tcd.1kg = SOURCE.LEVELS[c(2, 5)],
  all.datatypes.adx.asw = SOURCE.LEVELS
  )
KINSHIP.BIN.WIDTH <- 0.01
KINSHIP.DOWNSAMPLE.SIZES <- c(AFR = 118L, ADX = 50L, EUR = 119L)
BOOTSTRAP.REPLICATES <- 1000L
RANDOM.SEED <- 123L
PLOT.BASE.SIZE <- 24
CATEGORICAL.BAR.DODGE <- 0.9
CATEGORICAL.BAR.WIDTH <- 0.8
CATEGORICAL.BAR.LINEWIDTH <- 1
DENSE.BAR.WIDTH.MULTIPLIER <- 0.8
DENSE.BAR.LINEWIDTH <- 0.75
PLOT.STYLES <- list(
  population.colors = c(
    AFR = "#56B4E9", ADX = "#4B1FA8", EUR = "#fb8072",
    YRI = "#eec4dc", ASW = "#e44b8d", CEU = "#bb437e"
    ),
  source.colors = c(
    Simulation_2T12Consistent = "#9A83CE",
    Simulation_2T12Consistent_simDown = "#6F55B5",
    Simulation_largeGrowth = "#32146F",
    Simulation_largeGrowth_simDown = "#4B1FA8",
    Empirical = "#B83264"
    ),
  series.labels = SOURCE.LABELS
  )
BOOTSTRAP.PLOT.STYLES <- list(
  tcd.1kg = list(
    fill.colors = c(
      "T.C.D. AFR" = "#56B4E9", "T.C.D. ADX" = "#6F55B5",
      "T.C.D. EUR" = "#FB8072", YRI = "#EEC4DC", ASW = "#E44B8D",
      CEU = "#BB437E"
      ),
    fill.labels = c(
      "T.C.D. AFR" = "T.C.D. AFR", "T.C.D. ADX" = "T.C.D. ADX",
      "T.C.D. EUR" = "T.C.D. EUR", YRI = "YRI", ASW = "ASW",
      CEU = "CEU"
      )
    ),
  all.datatypes.adx.asw = list(
    fill.colors = c(
      "T.C." = "#9A83CE", "T.C.D." = "#6F55B5",
      "L.G." = "#32146F", "L.G.D." = "#4B1FA8", ASW = "#E44B8D"
      ),
    fill.labels = c(
      "T.C." = "T.C.", "T.C.D." = "T.C.D.", "L.G." = "L.G.",
      "L.G.D." = "L.G.D.", ASW = "ASW"
      )
    )
  )


# internal functions ----


# return the canonical levels represented by a filtered plot view
order.active.levels <- function(values, canonical.levels) {
  active.levels <- canonical.levels[
    canonical.levels %in% as.character(values)
    ]
  return(active.levels)
  }


# summarize the 50 complete simulation replicates with percentile intervals
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


# select population-specific empirical-sized simulation identifier sets
select.bootstrap.kinship.ids <- function(data, seed) {
  required.roles <- names(KINSHIP.DOWNSAMPLE.SIZES)
  id.column <- if ("sample.id" %in% names(data)) {
    "sample.id"
    } else if ("endpoint.1" %in% names(data)) {
    "endpoint.1"
    } else {
    stop("Kinship simulation data require an identifier column")
    }
  candidates <- data %>%
    filter(data.type != "Empirical", role %in% required.roles)
  candidates <- if (all(c("endpoint.1", "endpoint.2") %in% names(candidates))) {
    candidates %>%
      pivot_longer(c(endpoint.1, endpoint.2), values_to = "sample.id")
    } else {
    candidates %>% transmute(data.type, rep, chrom, role,
      sample.id = .data[[id.column]])
    }
  candidates <- candidates %>%
    select(data.type, rep, chrom, role, sample.id) %>%
    distinct()
  complete.role.sources <- candidates %>%
    filter(!grepl("largeGrowth", data.type)) %>%
    distinct(data.type, role) %>%
    count(data.type, name = "role.count")
  if (any(complete.role.sources$role.count != length(required.roles))) {
    stop("Kinship simulation data are missing required population roles")
    }
  sizes <- candidates %>% count(data.type, rep, chrom, role, name = "available")
  expected <- candidates %>% distinct(data.type, rep, chrom) %>%
    rowwise() %>% mutate(role = list(if (grepl("largeGrowth", data.type)) {
      "ADX" } else { required.roles })) %>% unnest(role)
  sample.sizes <- expected %>%
    left_join(sizes, by = c("data.type", "rep", "chrom", "role")) %>%
    mutate(available = replace_na(available, 0L),
      target = KINSHIP.DOWNSAMPLE.SIZES[role],
      selected = pmin(available, target),
      shortfall = target - selected)
  set.seed(seed)
  selected <- candidates %>%
    group_by(data.type, rep, chrom, role) %>%
    group_modify(function(group, key) {
      selected.size <- sample.sizes %>%
        filter(
          data.type == key$data.type,
          rep == key$rep,
          chrom == key$chrom,
          role == key$role
          ) %>%
        pull(selected)
      return(slice_sample(group, n = selected.size, replace = FALSE))
      }) %>%
    ungroup() %>%
    left_join(sample.sizes, by = c("data.type", "rep", "chrom", "role"))
  attr(selected, "sample.sizes") <- sample.sizes
  return(selected)
  }


# summarize per-replicate relationship fractions with simulation intervals
summarize.bootstrap.kinship <- function(data, breaks) {
  selected <- select.bootstrap.kinship.ids(data, RANDOM.SEED)
  selected.pairs <- apply.kinship.selection(data, selected) %>%
    filter(sample.set == "downsampled")
  histograms <- build.kinship.histograms(selected.pairs, breaks) %>%
    left_join(
      attr(selected, "sample.sizes"),
      by = c("data.type", "rep", "chrom", "role")
      )
  simulation <- histograms %>%
    filter(data.type != "Empirical") %>%
    summarize.bootstrap.interval(
      c(
        "data.type", "rep", "pop", "role", "chrom", "xmin", "xmax", "xmid"
        ),
      "fraction"
      )
  empirical <- summarize.empirical.kinship.interval(
    filter(data, data.type == "Empirical"), breaks,
    BOOTSTRAP.REPLICATES, RANDOM.SEED
    )
  summary <- bind_rows(simulation, empirical)
  attr(summary, "sample.sizes") <- attr(selected, "sample.sizes")
  return(summary)
  }


# resample empirical kinship pairs to calculate fixed-bin percentile intervals
summarize.empirical.kinship.interval <- function(
    data, breaks, replicates, seed
  ) {
  if (!is.numeric(replicates) || length(replicates) != 1L ||
      !is.finite(replicates) || replicates < 1L || replicates %% 1L != 0) {
    stop("Bootstrap replicate count must be a positive integer")
    }
  if (any(!is.finite(data$kinship))) {
    stop("Empirical kinship values must be finite")
    }
  set.seed(seed)
  summary <- data %>%
    group_by(data.type, pop, role, chrom) %>%
    group_modify(function(group, key) {
      observed <- hist(
        group$kinship, breaks = breaks, plot = FALSE, include.lowest = TRUE
        )
      observed.fraction <- observed$counts / sum(observed$counts)
      bootstrap.fractions <- vapply(
        seq_len(replicates),
        function(replicate.id) {
          values <- sample(group$kinship, nrow(group), replace = TRUE)
          bootstrap <- hist(
            values, breaks = breaks, plot = FALSE, include.lowest = TRUE
            )
          return(bootstrap$counts / sum(bootstrap$counts))
          },
        numeric(length(observed$counts))
        )
      bootstrap.fractions <- matrix(
        bootstrap.fractions, nrow = length(observed$counts)
        )
      return(tibble(
        xmin = head(observed$breaks, -1),
        xmax = tail(observed$breaks, -1),
        xmid = observed$mids,
        mean = observed.fraction,
        lower = apply(
          bootstrap.fractions, 1L, quantile, 0.025, names = FALSE
          ),
        upper = apply(
          bootstrap.fractions, 1L, quantile, 0.975, names = FALSE
          ),
        replicate.count = as.integer(replicates)
        ))
      }) %>%
    ungroup()
  return(summary)
  }


# select the source and population series for one bootstrap kinship view
filter.bootstrap.kinship.plot.view <- function(data, data.types, tag) {
  if (!tag %in% names(BOOTSTRAP.PLOT.STYLES)) {
    stop("Unsupported bootstrap kinship plot tag: ", tag)
    }
  if (!identical(data.types, PLOT.CONFIGS[[tag]])) {
    stop("Kinship data types do not match the configured bootstrap view")
    }
  plotted <- data %>%
    filter(
      as.character(data.type) %in% data.types,
      as.character(chrom) %in% SELECTED.CHROMOSOMES
      )
  plotted <- if (tag == "tcd.1kg") {
    plotted %>% filter(
      (data.type == "Simulation_2T12Consistent_simDown" &
         pop %in% c("AFR", "ADX", "EUR")) |
        (data.type == "Empirical" & pop %in% c("YRI", "ASW", "CEU"))
      )
    } else {
    plotted %>% filter(
      (data.type != "Empirical" & pop == "ADX") |
        (data.type == "Empirical" & pop == "ASW")
      )
    }
  styles <- BOOTSTRAP.PLOT.STYLES[[tag]]
  plotted <- plotted %>%
    mutate(
      data.type = factor(
        as.character(data.type),
        levels = order.active.levels(data.type, SOURCE.LEVELS)
        ),
      pop = factor(
        as.character(pop),
        levels = order.active.levels(pop, POPULATION.LEVELS)
        ),
      plot.key = case_when(
        tag == "tcd.1kg" & data.type != "Empirical" ~
          paste("T.C.D.", pop),
        tag == "tcd.1kg" ~ as.character(pop),
        data.type == "Empirical" ~ "ASW",
        TRUE ~ SOURCE.LABELS[as.character(data.type)]
        ),
      plot.key = factor(
        plot.key,
        levels = order.active.levels(plot.key, names(styles$fill.colors))
        )
      )
  return(plotted)
  }


# construct TCD/1kG and all-datatype ADX/ASW kinship plots with intervals
make.bootstrap.kinship.plot <- function(
    data, breaks, data.types, title, x.limits, tag
  ) {
  if (length(x.limits) != 2L || any(!is.finite(x.limits))) {
    stop("Kinship x limits must contain two finite values")
    }
  plotted <- filter.bootstrap.kinship.plot.view(data, data.types, tag)
  styles <- BOOTSTRAP.PLOT.STYLES[[tag]]
  dodge <- position_dodge(diff(breaks)[1])
  plot <- ggplot(plotted, aes(
    xmid, mean, fill = plot.key, group = plot.key
    )) +
    geom_col(
      position = dodge,
      width = diff(breaks)[1] * DENSE.BAR.WIDTH.MULTIPLIER,
      color = "black", linewidth = DENSE.BAR.LINEWIDTH
      ) +
    geom_errorbar(aes(ymin = lower, ymax = upper),
      position = dodge, width = 0, linewidth = DENSE.BAR.LINEWIDTH,
      na.rm = TRUE) +
    coord_cartesian(xlim = x.limits) +
    scale_fill_manual(
      values = styles$fill.colors,
      breaks = levels(plotted$plot.key),
      labels = styles$fill.labels[levels(plotted$plot.key)]
      ) +
    labs(title = title, x = "Pairwise KING kinship", y = "Fraction of pairs",
      fill = NULL) +
    guides(fill = guide_legend(order = 1, nrow = 1, byrow = TRUE)) +
    theme_bw(base_size = PLOT.BASE.SIZE) +
    theme(legend.position = "top", panel.grid.minor = element_blank(),)
  return(plot)
  }


# retain source-specific populations before any histogram calculations
apply.kinship.source.contract <- function(data) {
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
    stop("Kinship data contain an unsupported population label")
    }

  return(data)
  }


# standardize pair endpoints and metadata for one kinship table
normalize.kinship.table <- function(data, data.type.input) {
  if (all(c("sample1", "sample2") %in% names(data))) {
    endpoint.columns <- c("sample1", "sample2")
    } else if (all(c("id1", "id2") %in% names(data))) {
    endpoint.columns <- c("id1", "id2")
    } else {
    stop("Kinship data require sample1/sample2 or id1/id2 columns")
    }
  data <- data %>%
    mutate(
      rep = as.numeric(rep),
      chrom = as.character(chrom),
      pop = as.character(pop),
      endpoint.1 = as.character(.data[[endpoint.columns[1]]]),
      endpoint.2 = as.character(.data[[endpoint.columns[2]]]),
      kinship = as.numeric(kinship),
      data.type = data.type.input
      ) %>%
    select(
      rep, chrom, pop, endpoint.1, endpoint.2, kinship, data.type
      ) %>%
    apply.kinship.source.contract() %>%
    add.population.roles()

  return(data)
  }


# read chromosome-labelled kinship files for one source
read.kinship.chromosomes <- function(
    data.directory, chromosomes, data.type.input
  ) {
  paths <- file.path(
    path.expand(data.directory),
    glue::glue("kinship_unrelated.chr{chromosomes}.parquet")
    )
  missing <- !file.exists(paths)
  if (any(missing)) {
    warning(
      paste0(
        data.type.input,
        " kinship files are unavailable for chromosomes: ",
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
    return(normalize.kinship.table(table, data.type.input))
    })

  return(data)
  }


# select a reproducible fixed number of simulation identifiers per group
select.kinship.ids <- function(data, downsample.size, seed) {
  candidates <- data %>%
    filter(data.type != "Empirical") %>%
    select(data.type, rep, role, chrom, endpoint.1, endpoint.2) %>%
    pivot_longer(
      c(endpoint.1, endpoint.2), values_to = "sample.id"
      ) %>%
    distinct(data.type, rep, role, chrom, sample.id)
  expected.groups <- data %>%
    filter(data.type != "Empirical") %>%
    distinct(data.type, rep, chrom) %>%
    crossing(role = c("AFR", "ADX", "EUR"))
  sizes <- candidates %>%
    count(data.type, rep, role, chrom, name = "available")
  sizes <- expected.groups %>%
    left_join(
      sizes, by = c("data.type", "rep", "role", "chrom")
      ) %>%
    mutate(available = replace_na(available, 0L))
  if (any(sizes$available < downsample.size)) {
    stop(glue(
      "A simulation kinship group contains fewer than ",
      "{downsample.size} unique IDs"
      ))
    }
  set.seed(seed)
  selected <- candidates %>%
    group_by(data.type, rep, role, chrom) %>%
    slice_sample(n = downsample.size, replace = FALSE) %>%
    ungroup() %>%
    arrange(data.type, rep, role, chrom, sample.id)

  return(selected)
  }


# retain pairs only when both endpoints are in the selected identifier set
apply.kinship.selection <- function(data, selected.ids) {
  keys <- c("data.type", "rep", "role", "chrom")
  selected.endpoint.1 <- selected.ids %>%
    rename(endpoint.1 = sample.id)
  selected.endpoint.2 <- selected.ids %>%
    rename(endpoint.2 = sample.id)
  full <- data %>% mutate(sample.set = "full")
  downsampled <- data %>%
    inner_join(selected.endpoint.1, by = c(keys, "endpoint.1")) %>%
    inner_join(selected.endpoint.2, by = c(keys, "endpoint.2")) %>%
    mutate(sample.set = "downsampled")
  expected.counts <- selected.ids %>%
    count(across(all_of(keys)), name = "selected.count") %>%
    mutate(expected.pairs = choose(selected.count, 2))
  canonical.pairs <- downsampled %>%
    mutate(
      endpoint.low = pmin(endpoint.1, endpoint.2),
      endpoint.high = pmax(endpoint.1, endpoint.2)
      )
  actual.counts <- canonical.pairs %>%
    count(across(all_of(keys)), name = "actual.rows")
  unique.counts <- canonical.pairs %>%
    distinct(
      across(all_of(keys)), endpoint.low, endpoint.high
      ) %>%
    count(across(all_of(keys)), name = "actual.pairs")
  incomplete <- expected.counts %>%
    left_join(actual.counts, by = keys) %>%
    left_join(unique.counts, by = keys) %>%
    mutate(
      actual.rows = replace_na(actual.rows, 0L),
      actual.pairs = replace_na(actual.pairs, 0L)
      ) %>%
    filter(
      actual.rows != expected.pairs | actual.pairs != expected.pairs
      )
  has.self.pairs <- any(
    canonical.pairs$endpoint.1 == canonical.pairs$endpoint.2
    )
  if (nrow(incomplete) || has.self.pairs) {
    stop("Downsampled kinship groups require a complete pair table")
    }
  combined <- bind_rows(full, downsampled)

  return(combined)
  }


# define common histogram breaks spanning all finite kinship values
make.kinship.breaks <- function(data, bin.width) {
  values <- data$kinship[is.finite(data$kinship)]
  if (!length(values)) stop("Kinship data contain no finite values")
  lower <- floor(min(values) / bin.width) * bin.width
  upper <- ceiling(max(values) / bin.width) * bin.width
  if (lower == upper) upper <- lower + bin.width
  breaks <- seq(lower, upper, by = bin.width)

  return(breaks)
  }


# calculate source-normalized histogram fractions for every analysis group
build.kinship.histograms <- function(data, breaks) {
  histograms <- data %>%
    group_by(data.type, rep, pop, role, chrom) %>%
    group_modify(function(group, key) {
      histogram <- hist(
        group$kinship, breaks = breaks, plot = FALSE,
        include.lowest = TRUE
        )
      bins <- tibble(
        xmin = head(histogram$breaks, -1),
        xmax = tail(histogram$breaks, -1),
        xmid = histogram$mids,
        count = histogram$counts,
        fraction = histogram$counts / sum(histogram$counts),
        sample.size = n_distinct(c(group$endpoint.1, group$endpoint.2))
        )
      return(bins)
      }) %>%
    ungroup()

  return(histograms)
  }


# summarize simulation replicate histograms and retain empirical estimates
summarize.kinship.histograms <- function(data) {
  simulation <- data %>%
    filter(data.type != "Empirical") %>%
    group_by(
      data.type, pop, role, chrom, xmin, xmax, xmid
      ) %>%
    summarise(
      mean.fraction = mean(fraction),
      sd.fraction = sd(fraction),
      replicate.count = n_distinct(rep),
      sample.size.min = min(sample.size),
      sample.size.max = max(sample.size),
      .groups = "drop"
      )
  empirical <- data %>%
    filter(data.type == "Empirical") %>%
    transmute(
      data.type, pop, role, chrom, xmin, xmax, xmid,
      mean.fraction = fraction,
      sample.size.min = sample.size,
      sample.size.max = sample.size
      ) %>%
    distinct() %>%
    mutate(sd.fraction = NA_real_, replicate.count = 1L)
  summary <- bind_rows(simulation, empirical) %>%
    mutate(
      pop = factor(pop, levels = POPULATION.LEVELS),
      chrom = factor(
        chrom, levels = c(SELECTED.CHROMOSOMES, "all")
        ),
      data.type = factor(
        data.type,
        levels = SOURCE.LEVELS
        )
      )

  return(summary)
  }


# retain one configured pairwise kinship view
filter.plot.view <- function(data, data.types, tag) {
  if (!tag %in% names(PLOT.CONFIGS)) {
    stop("Unsupported kinship plot tag: ", tag)
    }
  if (!identical(data.types, PLOT.CONFIGS[[tag]])) {
    stop("Kinship data types do not match the configured tag")
    }
  filtered <- data %>%
    filter(as.character(data.type) %in% data.types) %>%
    filter(tag != "onlyADX" | pop == "ADX")
  active.sources <- order.active.levels(filtered$data.type, SOURCE.LEVELS)
  active.populations <- order.active.levels(
    filtered$pop, POPULATION.LEVELS
    )
  filtered <- filtered %>%
    mutate(
      data.type = factor(as.character(data.type), levels = active.sources),
      pop = factor(as.character(pop), levels = active.populations)
      ) %>%
    arrange(data.type) %>%
    droplevels()
  return(filtered)
  }


# construct one scoped pairwise kinship distribution plot
make.kinship.plot <- function(
    data, breaks, styles, data.types, tag, show.all = TRUE
  ) {
  source.view <- tag == "onlyADX"
  plot.data <- filter.plot.view(data, data.types, tag) %>%
    filter(
      as.character(chrom) %in% SELECTED.CHROMOSOMES |
        (show.all & data.type == "Empirical" & chrom == "all")
      ) %>%
    mutate(chrom = factor(
      as.character(chrom),
      levels = c(SELECTED.CHROMOSOMES, if (show.all) "all")
      )) %>%
    mutate(
      plot.key = if (source.view) {
        factor(
          as.character(data.type),
          levels = order.active.levels(data.type, SOURCE.LEVELS)
          )
        } else {
        factor(
          as.character(pop),
          levels = order.active.levels(pop, POPULATION.LEVELS)
          )
        }
      )
  dodge <- position_dodge(width = diff(breaks)[1] * 0.9)
  plot <- ggplot(
    plot.data,
    aes(
      x = xmid, y = mean.fraction, fill = plot.key,
      group = plot.key
      )
    ) +
    geom_col(
      position = dodge, width = diff(breaks)[1] * 0.85,
      color = "black", linewidth = 0.1
      ) +
    geom_errorbar(
      data = plot.data,
      aes(
        ymin = pmax(0, mean.fraction - 2 * sd.fraction),
        ymax = mean.fraction + 2 * sd.fraction
        ),
      position = dodge, width = 0, linewidth = 0.45,
      color = "black", na.rm = TRUE
      ) +
    facet_grid(
      chrom ~ data.type, drop = TRUE, scales = "free_y",
      labeller = labeller(
        data.type = as_labeller(styles$series.labels)
        )
      ) +
    scale_fill_manual(
      values = if (source.view) styles$source.colors else {
        styles$population.colors
        },
      breaks = levels(plot.data$plot.key),
      labels = if (source.view) {
        styles$series.labels[levels(plot.data$plot.key)]
        } else {
        waiver()
        }
      ) +
    labs(
      x = "Pairwise KING Kinship", y = "Fraction of pairs",
      title = paste("Pairwise KING Kinship Distributions:", tag),
      fill = NULL
      ) +
    xlim(-0.2, 0.0442) +
    guides(fill = guide_legend(order = 1, nrow = 1, byrow = TRUE)) +
    theme_bw(base_size = PLOT.BASE.SIZE) +
    theme(
      legend.position = "top",
      legend.direction = "horizontal",
      legend.box = "horizontal",
      legend.title = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold")
      )

  return(plot)
  }


# analysis ----


# read selected chromosome-level simulation and empirical estimates
sim.tc.kinship <- read.kinship.chromosomes(
  SIM.TC.DATA.DIR, SELECTED.CHROMOSOMES, "Simulation_2T12Consistent"
  )
simDown.tc.kinship <- read.kinship.chromosomes(
  SIMDOWN.TC.DATA.DIR, SELECTED.CHROMOSOMES,
  "Simulation_2T12Consistent_simDown"
  )
sim.lg.kinship <- read.kinship.chromosomes(
  SIM.LG.DATA.DIR, SELECTED.CHROMOSOMES, "Simulation_largeGrowth"
  )
simDown.lg.kinship <- read.kinship.chromosomes(
  SIMDOWN.LG.DATA.DIR, SELECTED.CHROMOSOMES,
  "Simulation_largeGrowth_simDown"
  )
emp.chromosome.kinship <- read.kinship.chromosomes(
  EMPIRICAL.DATA.DIR, SELECTED.CHROMOSOMES, "Empirical"
  )

# combine all unrelated simulation and empirical pairs
simulation.kinship <- bind_rows(
  sim.tc.kinship, simDown.tc.kinship,
  sim.lg.kinship, simDown.lg.kinship
  )
empirical.kinship <- emp.chromosome.kinship

# summarize common-bin histograms and construct all configured plots
kinship.data <- bind_rows(simulation.kinship, empirical.kinship)
kinship.breaks <- make.kinship.breaks(
  kinship.data, KINSHIP.BIN.WIDTH
  )
kinship.summary <- summarize.bootstrap.kinship(kinship.data, kinship.breaks)
kinship.bootstrap.tcd.1kg <- make.bootstrap.kinship.plot(
  kinship.summary, kinship.breaks, PLOT.CONFIGS$tcd.1kg,
  "Pairwise KING kinship: chromosome 1 TCD and 1kG",
  x.limits = c(-0.1, 0.05), tag = "tcd.1kg"
  )
kinship.bootstrap.all.datatypes.adx.asw <- make.bootstrap.kinship.plot(
  kinship.summary, kinship.breaks, PLOT.CONFIGS$all.datatypes.adx.asw,
  "Pairwise KING kinship: chromosome 1 all ADX sources and ASW",
  x.limits = c(-0.2, 0.05), tag = "all.datatypes.adx.asw"
  )

# persist bootstrap plots before printing figures at the end of the script
dir.create(OUTPUT.DIR, recursive = TRUE, showWarnings = FALSE)
write_csv(attr(kinship.summary, "sample.sizes"), file.path(
  OUTPUT.DIR, "kinship.bootstrap.sample_sizes.csv"
  ))
saveRDS(kinship.bootstrap.tcd.1kg, file.path(
  OUTPUT.DIR, "kinship.bootstrap.tcd.1kg.rds"
  ))
saveRDS(kinship.bootstrap.all.datatypes.adx.asw, file.path(
  OUTPUT.DIR, "kinship.bootstrap.all.datatypes.adx.asw.rds"
  ))

print(kinship.bootstrap.tcd.1kg)
print(kinship.bootstrap.all.datatypes.adx.asw)

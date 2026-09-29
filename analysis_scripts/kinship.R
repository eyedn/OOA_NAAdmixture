# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology
# Mooney Lab
# ___
# kinship.R
# ______________________________________________________________________________


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
  tcd.1kg = SOURCE.LEVELS[c(2, 5)],
  tc.tcd.1kg = SOURCE.LEVELS[c(1, 2, 5)],
  all.datatypes.adx.asw = SOURCE.LEVELS
  )
KINSHIP.BIN.WIDTH <- 0.01
KINSHIP.DOWNSAMPLE.SIZES <- c(AFR = 118L, ADX = 50L, EUR = 119L)
BOOTSTRAP.REPLICATES <- 1000L
RANDOM.SEED <- 123L
PLOT.BASE.SIZE <- 24
DENSE.BAR.WIDTH.MULTIPLIER <- 0.8
DENSE.BAR.LINEWIDTH <- 0.75
KINSHIP.CONTRAST.BOOTSTRAP.REPLICATES <- 50L
KINSHIP.CONTRAST.FAMILY.SIZE <- 75L
KINSHIP.CONTRAST.X.LIMITS <- c(-0.20, 0.05)
KINSHIP.POPULATION.CONTRAST.SOURCES <- SOURCE.LEVELS[1:2]
KINSHIP.POPULATION.CONTRASTS <- tribble(
  ~contrast, ~left.pop, ~right.pop,
  "AFR-ADX", "AFR", "ADX",
  "AFR-EUR", "AFR", "EUR",
  "ADX-EUR", "ADX", "EUR"
  )
KINSHIP.SIMULATION.CONTRAST.SOURCES <- SOURCE.LEVELS[1:4]
KINSHIP.SIMULATION.CONTRASTS <- tribble(
  ~contrast, ~left.source, ~right.source,
  "T.C. - T.C.D.", SOURCE.LEVELS[[1L]], SOURCE.LEVELS[[2L]],
  "L.G. - L.G.D.", SOURCE.LEVELS[[3L]], SOURCE.LEVELS[[4L]],
  "T.C. - L.G.", SOURCE.LEVELS[[1L]], SOURCE.LEVELS[[3L]]
  )
KINSHIP.EMPIRICAL.CONTRASTS <- tribble(
  ~contrast, ~left.pop, ~right.pop,
  "YRI-ASW", "YRI", "ASW",
  "YRI-CEU", "YRI", "CEU",
  "ASW-CEU", "ASW", "CEU"
  )
KINSHIP.CONTRAST.COLORS <- c(
  "AFR-ADX" = "#BDBDBD", "AFR-EUR" = "#737373", "ADX-EUR" = "#000000",
  "T.C. - T.C.D." = "#BDBDBD", "L.G. - L.G.D." = "#737373",
  "T.C. - L.G." = "#000000",
  "YRI-ASW" = "#BDBDBD", "YRI-CEU" = "#737373", "ASW-CEU" = "#000000"
  )
KINSHIP.CONTRAST.LABELS <- c(
  "AFR-ADX" = "AFR - ADX", "AFR-EUR" = "AFR - EUR",
  "ADX-EUR" = "ADX - EUR", "T.C. - T.C.D." = "T.C. - T.C.D.",
  "L.G. - L.G.D." = "L.G. - L.G.D.", "T.C. - L.G." = "T.C. - L.G.",
  "YRI-ASW" = "YRI - ASW", "YRI-CEU" = "YRI - CEU",
  "ASW-CEU" = "ASW - CEU"
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


# select exactly the displayed bins while retaining all-pair normalization
select.kinship.contrast.bins <- function(data) {
  selected <- data %>% filter(
    xmin >= KINSHIP.CONTRAST.X.LIMITS[[1L]] - .Machine$double.eps,
    xmax <= KINSHIP.CONTRAST.X.LIMITS[[2L]] + .Machine$double.eps
    )
  expected.bins <- diff(KINSHIP.CONTRAST.X.LIMITS) / KINSHIP.BIN.WIDTH
  if (n_distinct(selected$xmin) != expected.bins) {
    stop("Kinship contrast display requires exactly 25 fixed bins")
    }
  return(selected)
  }


# require complete finite per-replicate simulation fractions for contrasts
validate.kinship.contrast.input <- function(data, sources, populations) {
  required.columns <- c("data.type", "rep", "chrom", "pop", "xmin",
    "xmax", "xmid", "fraction")
  missing.columns <- setdiff(required.columns, names(data))
  if (length(missing.columns)) {
    stop("Kinship contrast inputs are missing columns: ",
      paste(missing.columns, collapse = ", "))
    }
  data <- data %>%
    filter(data.type %in% sources, pop %in% populations,
      as.character(chrom) %in% SELECTED.CHROMOSOMES) %>%
    select.kinship.contrast.bins() %>%
    mutate(data.type = as.character(data.type), chrom = as.character(chrom),
      pop = as.character(pop))
  expected <- crossing(data.type = sources, rep = seq_len(50L),
    chrom = SELECTED.CHROMOSOMES, pop = populations,
    xmin = sort(unique(data$xmin)))
  counts <- data %>% count(data.type, rep, chrom, pop, xmin,
    name = "fraction.count")
  invalid <- expected %>% left_join(counts,
    by = c("data.type", "rep", "chrom", "pop", "xmin")) %>%
    mutate(fraction.count = replace_na(fraction.count, 0L)) %>%
    filter(fraction.count != 1L)
  if (nrow(invalid) || any(!is.finite(data$fraction))) {
    stop("Kinship contrast inputs require exactly one finite fraction for ",
      "every source, replicate, population, and bin")
    }
  return(data)
  }


# bootstrap a paired left-minus-right kinship fraction difference
summarize.kinship.contrast.bootstrap <- function(
    left.values, right.values, bootstrap.replicates,
    family.size = KINSHIP.CONTRAST.FAMILY.SIZE
  ) {
  if (length(left.values) != 50L || length(right.values) != 50L ||
      any(!is.finite(left.values)) || any(!is.finite(right.values))) {
    stop("Kinship contrast bootstraps require 50 finite paired values")
    }
  draws <- replicate(bootstrap.replicates, {
    indices <- sample(seq_along(left.values), length(left.values),
      replace = TRUE)
    mean(left.values[indices] - right.values[indices])
    })
  nominal <- quantile(draws, c(0.025, 0.975), names = FALSE)
  bonferroni.quantile <- 0.05 / (2 * family.size)
  bonferroni <- quantile(draws,
    c(bonferroni.quantile, 1 - bonferroni.quantile), names = FALSE)
  return(tibble(
    difference = mean(left.values) - mean(right.values),
    ci.95.lower = nominal[[1L]], ci.95.upper = nominal[[2L]],
    bonferroni.ci.lower = bonferroni[[1L]],
    bonferroni.ci.upper = bonferroni[[2L]],
    bonferroni.quantile = bonferroni.quantile
    ))
  }


# construct paired population contrast intervals from simulation histograms
make.kinship.population.contrast.tables <- function(
    data, bootstrap.replicates = KINSHIP.CONTRAST.BOOTSTRAP.REPLICATES,
    seed = RANDOM.SEED
  ) {
  data <- validate.kinship.contrast.input(data,
    KINSHIP.POPULATION.CONTRAST.SOURCES, c("AFR", "ADX", "EUR"))
  set.seed(seed)
  tables <- map_dfr(KINSHIP.POPULATION.CONTRAST.SOURCES, function(source) {
    map_dfr(seq_len(nrow(KINSHIP.POPULATION.CONTRASTS)), function(index) {
      contrast <- KINSHIP.POPULATION.CONTRASTS[index, ]
      map_dfr(sort(unique(data$xmin)), function(bin) {
        left <- data %>% filter(data.type == source,
          pop == contrast$left.pop, xmin == bin) %>% arrange(rep)
        right <- data %>% filter(data.type == source,
          pop == contrast$right.pop, xmin == bin) %>% arrange(rep)
        if (!identical(left$rep, right$rep)) {
          stop("Paired kinship contrast replicate IDs must match")
          }
        bind_cols(tibble(data.type = source, xmin = bin,
          xmax = left$xmax[[1L]], xmid = left$xmid[[1L]],
          contrast = contrast$contrast),
          summarize.kinship.contrast.bootstrap(left$fraction, right$fraction,
            bootstrap.replicates))
        })
      })
    })
  return(tables)
  }


# construct paired ADX simulation-source contrast intervals from histograms
make.kinship.simulation.contrast.tables <- function(
    data, bootstrap.replicates = KINSHIP.CONTRAST.BOOTSTRAP.REPLICATES,
    seed = RANDOM.SEED
  ) {
  data <- validate.kinship.contrast.input(data,
    KINSHIP.SIMULATION.CONTRAST.SOURCES, "ADX")
  set.seed(seed)
  tables <- map_dfr(
    seq_len(nrow(KINSHIP.SIMULATION.CONTRASTS)), function(index) {
    contrast <- KINSHIP.SIMULATION.CONTRASTS[index, ]
    map_dfr(sort(unique(data$xmin)), function(bin) {
      left <- data %>% filter(data.type == contrast$left.source,
        xmin == bin) %>% arrange(rep)
      right <- data %>% filter(data.type == contrast$right.source,
        xmin == bin) %>% arrange(rep)
      if (!identical(left$rep, right$rep)) {
        stop("Paired kinship contrast replicate IDs must match")
        }
      bind_cols(tibble(xmin = bin, xmax = left$xmax[[1L]],
        xmid = left$xmid[[1L]], contrast = contrast$contrast),
        summarize.kinship.contrast.bootstrap(left$fraction, right$fraction,
          bootstrap.replicates))
      })
      })
  return(tables)
  }


# bootstrap independently resampled empirical pair fractions for each contrast
make.kinship.empirical.contrast.tables <- function(
    data, breaks, bootstrap.replicates = KINSHIP.CONTRAST.BOOTSTRAP.REPLICATES,
    seed = RANDOM.SEED
  ) {
  if (any(!is.finite(data$kinship))) {
    stop("Empirical kinship contrast values must be finite")
    }
  data <- data %>% filter(data.type == "Empirical",
    as.character(chrom) %in% SELECTED.CHROMOSOMES,
    pop %in% c("YRI", "ASW", "CEU"))
  if (!all(c("YRI", "ASW", "CEU") %in% unique(data$pop))) {
    stop("Empirical kinship contrast inputs require YRI, ASW, and CEU pairs")
    }
  set.seed(seed)
  displayed.breaks <- breaks[
    breaks >= KINSHIP.CONTRAST.X.LIMITS[[1L]] - .Machine$double.eps &
      breaks <= KINSHIP.CONTRAST.X.LIMITS[[2L]] + .Machine$double.eps
    ]
  if (length(displayed.breaks) != 26L) {
    stop("Kinship contrast display requires exactly 25 fixed bins")
    }
  tables <- map_dfr(
    seq_len(nrow(KINSHIP.EMPIRICAL.CONTRASTS)), function(index) {
    contrast <- KINSHIP.EMPIRICAL.CONTRASTS[index, ]
    left <- filter(data, pop == contrast$left.pop)$kinship
    right <- filter(data, pop == contrast$right.pop)$kinship
    observed <- function(values) {
      histogram <- hist(values, breaks = breaks, plot = FALSE,
        include.lowest = TRUE)
      histogram$counts / length(values)
      }
    left.observed <- observed(left)
    right.observed <- observed(right)
    draws <- replicate(bootstrap.replicates, {
      observed(sample(left, length(left), replace = TRUE)) -
        observed(sample(right, length(right), replace = TRUE))
      })
    bin.index <- match(head(displayed.breaks, -1L), head(breaks, -1L))
    nominal <- apply(draws[bin.index, , drop = FALSE], 1L, quantile,
      c(0.025, 0.975), names = FALSE)
    bonferroni.quantile <- 0.05 / (2 * KINSHIP.CONTRAST.FAMILY.SIZE)
    bonferroni <- apply(draws[bin.index, , drop = FALSE], 1L, quantile,
      c(bonferroni.quantile, 1 - bonferroni.quantile), names = FALSE)
    tibble(
      xmin = head(displayed.breaks, -1L),
      xmax = tail(displayed.breaks, -1L),
      xmid = rowMeans(embed(displayed.breaks, 2L)),
      contrast = contrast$contrast,
      difference = left.observed[bin.index] - right.observed[bin.index],
      ci.95.lower = nominal[1L, ], ci.95.upper = nominal[2L, ],
      bonferroni.ci.lower = bonferroni[1L, ],
      bonferroni.ci.upper = bonferroni[2L, ],
      bonferroni.quantile = bonferroni.quantile)
      })
  return(tables)
  }


# select one interval family and preserve the configured contrast order
prepare.kinship.contrast.plot.data <- function(
    data, interval.type = c("95", "bonferroni"),
    family = c("simulation", "empirical")
  ) {
  interval.type <- match.arg(interval.type)
  family <- match.arg(family)
  lower <- if (interval.type == "95") "ci.95.lower" else "bonferroni.ci.lower"
  upper <- if (interval.type == "95") "ci.95.upper" else "bonferroni.ci.upper"
  contrasts <- if (family == "empirical") {
    KINSHIP.EMPIRICAL.CONTRASTS$contrast
    } else if ("data.type" %in% names(data)) {
    KINSHIP.POPULATION.CONTRASTS$contrast
    } else {
    KINSHIP.SIMULATION.CONTRASTS$contrast
    }
  displayed <- data %>% transmute(across(any_of("data.type")), xmin, xmax,
    xmid, contrast, difference, ci.lower = .data[[lower]],
    ci.upper = .data[[upper]]) %>% mutate(
      contrast = factor(contrast, levels = contrasts),
      outside.zero = ci.lower > 0 | ci.upper < 0,
      point.color = if_else(
        outside.zero, "red", as.character(contrast)
        )
      )
  if ("data.type" %in% names(displayed)) {
    displayed <- displayed %>% mutate(data.type = factor(data.type,
      levels = KINSHIP.POPULATION.CONTRAST.SOURCES))
    }
  return(list(
    simulation = displayed,
    red.markers = filter(displayed, outside.zero)
    ))
  }


# build a fixed-bin kinship contrast plot with zero reference and intervals
make.kinship.contrast.plot <- function(data, interval.label, title) {
  simulation <- data$simulation
  contrasts <- levels(simulation$contrast)
  dodge <- position_dodge(KINSHIP.BIN.WIDTH)
  plot <- ggplot(
    simulation, aes(xmid, difference, group = contrast)
    ) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_errorbar(
      aes(ymin = ci.lower, ymax = ci.upper, color = contrast),
      position = dodge, width = 0, linewidth = DENSE.BAR.LINEWIDTH,
      show.legend = FALSE
      ) +
    geom_point(
      aes(
        xmid, difference, color = point.color, fill = contrast,
        group = contrast
        ),
      shape = 21, position = dodge, size = 2.5,
      stroke = DENSE.BAR.LINEWIDTH, inherit.aes = FALSE
      ) +
    coord_cartesian(xlim = KINSHIP.CONTRAST.X.LIMITS) +
    scale_color_manual(
      values = c(KINSHIP.CONTRAST.COLORS, red = "red"), breaks = contrasts,
      labels = KINSHIP.CONTRAST.LABELS[contrasts]
      ) +
    scale_fill_manual(values = KINSHIP.CONTRAST.COLORS[contrasts],
      guide = "none") +
    labs(x = "Pairwise KING kinship", y = "Left minus right fraction of pairs",
      color = NULL, title = title,
      subtitle = paste(interval.label, "interval")) +
    guides(
      color = guide_legend(
        nrow = 1, byrow = TRUE,
        override.aes = list(
          shape = 21, fill = unname(KINSHIP.CONTRAST.COLORS[contrasts]),
          color = unname(KINSHIP.CONTRAST.COLORS[contrasts])
          )
        )
      ) +
    theme_bw(base_size = PLOT.BASE.SIZE) +
    theme(legend.position = "top", panel.grid.minor = element_blank())
  if ("data.type" %in% names(simulation)) {
    plot <- plot + facet_grid(
      . ~ data.type, labeller = labeller(data.type = SOURCE.LABELS)
      )
    }
  return(plot)
  }


# analysis ----


# read selected chromosome-level simulation and empirical estimates
kinship.sim.tc <- read.kinship.chromosomes(
  SIM.TC.DATA.DIR, SELECTED.CHROMOSOMES, "Simulation_2T12Consistent"
  )
kinship.simDown.tc <- read.kinship.chromosomes(
  SIMDOWN.TC.DATA.DIR, SELECTED.CHROMOSOMES,
  "Simulation_2T12Consistent_simDown"
  )
kinship.sim.lg <- read.kinship.chromosomes(
  SIM.LG.DATA.DIR, SELECTED.CHROMOSOMES, "Simulation_largeGrowth"
  )
kinship.simDown.lg <- read.kinship.chromosomes(
  SIMDOWN.LG.DATA.DIR, SELECTED.CHROMOSOMES,
  "Simulation_largeGrowth_simDown"
  )
kinship.emp.chromosome <- read.kinship.chromosomes(
  EMPIRICAL.DATA.DIR, SELECTED.CHROMOSOMES, "Empirical"
  )

# combine all unrelated simulation and empirical pairs
kinship.simulation <- bind_rows(
  kinship.sim.tc, kinship.simDown.tc,
  kinship.sim.lg, kinship.simDown.lg
  )
kinship.empirical <- kinship.emp.chromosome

# summarize common-bin histograms and construct all configured plots
kinship.data <- bind_rows(kinship.simulation, kinship.empirical)
kinship.breaks <- make.kinship.breaks(
  kinship.data, KINSHIP.BIN.WIDTH
  )
kinship.summary <- summarize.bootstrap.kinship(kinship.data, kinship.breaks)
kinship.contrast.selected.ids <- select.bootstrap.kinship.ids(
  kinship.data, RANDOM.SEED
  )
kinship.contrast.histograms <- apply.kinship.selection(
  kinship.data, kinship.contrast.selected.ids
  ) %>%
  filter(sample.set == "downsampled") %>%
  build.kinship.histograms(kinship.breaks)
kinship.population.contrast.tables <- make.kinship.population.contrast.tables(
  filter(kinship.contrast.histograms, data.type != "Empirical")
  )
kinship.simulation.contrast.tables <- make.kinship.simulation.contrast.tables(
  filter(kinship.contrast.histograms, data.type != "Empirical")
  )
kinship.empirical.contrast.tables <- make.kinship.empirical.contrast.tables(
  kinship.empirical, kinship.breaks
  )
kinship.population.contrast.plots <- list(
  `95` = make.kinship.contrast.plot(
    prepare.kinship.contrast.plot.data(
      kinship.population.contrast.tables, "95"
      ),
    "95%", "T.C. and T.C.D. kinship population contrasts"
    ),
  bonferroni = make.kinship.contrast.plot(
    prepare.kinship.contrast.plot.data(
      kinship.population.contrast.tables, "bonferroni"
      ),
    "Bonferroni", "T.C. and T.C.D. kinship population contrasts"
    )
  )
kinship.simulation.contrast.plots <- list(
  `95` = make.kinship.contrast.plot(
    prepare.kinship.contrast.plot.data(
      kinship.simulation.contrast.tables, "95"
      ),
    "95%", "ADX kinship simulation-source contrasts"
    ),
  bonferroni = make.kinship.contrast.plot(
    prepare.kinship.contrast.plot.data(
      kinship.simulation.contrast.tables, "bonferroni"
      ),
    "Bonferroni", "ADX kinship simulation-source contrasts"
    )
  )
kinship.empirical.contrast.plots <- list(
  `95` = make.kinship.contrast.plot(
    prepare.kinship.contrast.plot.data(
      kinship.empirical.contrast.tables, "95", "empirical"
      ),
    "95%", "Empirical kinship population contrasts"
    ),
  bonferroni = make.kinship.contrast.plot(
    prepare.kinship.contrast.plot.data(
      kinship.empirical.contrast.tables, "bonferroni", "empirical"
      ),
    "Bonferroni", "Empirical kinship population contrasts"
    )
  )
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
saveRDS(kinship.population.contrast.plots$`95`, file.path(
  OUTPUT.DIR, "kinship.population.contrasts.95.rds"
  ))
saveRDS(kinship.population.contrast.plots$bonferroni, file.path(
  OUTPUT.DIR, "kinship.population.contrasts.bonferroni.rds"
  ))
saveRDS(kinship.simulation.contrast.plots$`95`, file.path(
  OUTPUT.DIR, "kinship.simulation.contrasts.95.rds"
  ))
saveRDS(kinship.simulation.contrast.plots$bonferroni, file.path(
  OUTPUT.DIR, "kinship.simulation.contrasts.bonferroni.rds"
  ))
saveRDS(kinship.empirical.contrast.plots$`95`, file.path(
  OUTPUT.DIR, "kinship.empirical.contrasts.95.rds"
  ))
saveRDS(kinship.empirical.contrast.plots$bonferroni, file.path(
  OUTPUT.DIR, "kinship.empirical.contrasts.bonferroni.rds"
  ))

print(kinship.bootstrap.tcd.1kg)
print(kinship.bootstrap.all.datatypes.adx.asw)
print(kinship.population.contrast.plots$`95`)
print(kinship.population.contrast.plots$bonferroni)
print(kinship.simulation.contrast.plots$`95`)
print(kinship.simulation.contrast.plots$bonferroni)
print(kinship.empirical.contrast.plots$`95`)
print(kinship.empirical.contrast.plots$bonferroni)

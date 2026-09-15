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


SIM.TC.DATA.DIR <- "~/scratch/OOA_NAAdmixture_small/stats"
SIMDOWN.TC.DATA.DIR <- "~/scratch/OOA_NAAdmixture_smallOnekgDownsample/stats"
SIM.LG.DATA.DIR <- "~/scratch/OOA_NAAdmixture_large/stats"
SIMDOWN.LG.DATA.DIR <- "~/scratch/OOA_NAAdmixture_largeOnekgDownsample/stats"
EMPIRICAL.DATA.DIR <- "~/scratch/OOA_NAAdmixture_1kG/stats"
OUTPUT.DIR <- "/home1/karatas/proj/OOA_NAAdmixture_data"
CHROMOSOMES <- as.character(1:22)
SELECTED.CHROMOSOMES <- c("1", "18")
SOURCE.LEVELS <- c(
  "Simulation_2T12Consistent", "Simulation_2T12Consistent_simDown",
  "Simulation_largeGrowth", "Simulation_largeGrowth_simDown",
  "Empirical"
  )
PLOT.CONFIGS <- list(
  TC.1kG = SOURCE.LEVELS[c(1, 5)],
  TC.TCD = SOURCE.LEVELS[c(1, 2)],
  TCD.1kG = SOURCE.LEVELS[c(2, 5)],
  onlyADX = SOURCE.LEVELS[1:4]
  )
KINSHIP.BIN.WIDTH <- 0.01
PLOT.BASE.SIZE <- 24
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
  series.labels = c(
    Simulation_2T12Consistent = "TC",
    Simulation_2T12Consistent_simDown = "TC D.",
    Simulation_largeGrowth = "LG",
    Simulation_largeGrowth_simDown = "LG D.",
    Empirical = "Emp."
    )
  )


# internal functions ----


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


# read the empirical genome-wide kinship table
read.empirical.kinship.genome <- function(data.directory) {
  data <- read_parquet(file.path(
    path.expand(data.directory), "kinship_unrelated.parquet"
    ))
  data$chrom <- "all"
  data <- normalize.kinship.table(data, "Empirical")

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
      pop = factor(pop, levels = c("AFR", "ADX", "EUR", "YRI", "ASW", "CEU")),
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
    filter(tag != "onlyADX" | pop == "ADX") %>%
    mutate(data.type = factor(as.character(data.type), levels = data.types)) %>%
    arrange(data.type) %>%
    droplevels()
  return(filtered)
  }


# construct one scoped pairwise kinship distribution plot
make.kinship.plot <- function(
    data, breaks, styles, data.types, tag, show.all = FALSE
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
        factor(as.character(data.type), levels = data.types)
        } else {
        factor(as.character(pop), levels = c(
          "AFR", "ADX", "EUR", "YRI", "ASW", "CEU"
          ))
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
    scale_fill_manual(values = if (source.view) {
      styles$source.colors
      } else {
      styles$population.colors
      }, breaks = if (source.view) {
      data.types
      } else {
      names(styles$population.colors)
      }, labels = if (source.view) {
      styles$series.labels[data.types]
      } else {
      waiver()
      }) +
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


# read all available chromosome-level simulation and empirical estimates
sim.tc.kinship <- read.kinship.chromosomes(
  SIM.TC.DATA.DIR, CHROMOSOMES, "Simulation_2T12Consistent"
  )
simDown.tc.kinship <- read.kinship.chromosomes(
  SIMDOWN.TC.DATA.DIR, CHROMOSOMES,
  "Simulation_2T12Consistent_simDown"
  )
sim.lg.kinship <- read.kinship.chromosomes(
  SIM.LG.DATA.DIR, CHROMOSOMES, "Simulation_largeGrowth"
  )
simDown.lg.kinship <- read.kinship.chromosomes(
  SIMDOWN.LG.DATA.DIR, CHROMOSOMES,
  "Simulation_largeGrowth_simDown"
  )
emp.chromosome.kinship <- read.kinship.chromosomes(
  EMPIRICAL.DATA.DIR, CHROMOSOMES, "Empirical"
  )
emp.genome.kinship <- read.empirical.kinship.genome(
  EMPIRICAL.DATA.DIR
  )

# combine all unrelated simulation and empirical pairs
simulation.kinship <- bind_rows(
  sim.tc.kinship, simDown.tc.kinship,
  sim.lg.kinship, simDown.lg.kinship
  )
empirical.kinship <- bind_rows(
  emp.chromosome.kinship, emp.genome.kinship
  )

# summarize common-bin histograms and construct all configured plots
kinship.data <- bind_rows(simulation.kinship, empirical.kinship)
kinship.breaks <- make.kinship.breaks(
  kinship.data, KINSHIP.BIN.WIDTH
  )
kinship.histograms <- build.kinship.histograms(
  kinship.data, kinship.breaks
  )
kinship.summary <- summarize.kinship.histograms(
  kinship.histograms
  )
kinship.plots <- imap(PLOT.CONFIGS, function(data.types, tag) {
  return(make.kinship.plot(
    kinship.summary, kinship.breaks, PLOT.STYLES, data.types, tag,
    show.all = FALSE
    ))
  })

# persist every plot before printing figures at the end of the script
dir.create(OUTPUT.DIR, recursive = TRUE, showWarnings = FALSE)
iwalk(kinship.plots, function(plot, tag) {
  saveRDS(plot, file.path(
    OUTPUT.DIR,
    str_replace("kinship.{tag}.rds", fixed("{tag}"), tag)
    ))
  })

print(kinship.plots$TC.1kG)
print(kinship.plots$TC.TCD)
print(kinship.plots$TCD.1kG)
print(kinship.plots$onlyADX)

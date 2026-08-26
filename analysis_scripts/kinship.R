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
library(nanoparquet)


SIM.SMALL.DATA.DIR <- "~/scratch/OOA_NAAdmixture_small/stats"
SIM.LARGE.DATA.DIR <- "~/scratch/OOA_NAAdmixture_large/stats"
EMPIRICAL.DATA.DIR <- "~/scratch/OOA_NAAdmixture_1kG/stats"
CHROMOSOMES <- as.character(1:22)
SELECTED.CHROMOSOMES <- c("1", "18")
DOWNSAMPLE.SIZE <- 50
RANDOM.SEED <- 123
KINSHIP.BIN.WIDTH <- 0.01
PLOT.SAMPLE.SET <- "downsampled"
PLOT.BASE.SIZE <- 24
PLOT.STYLES <- list(
  population.colors = c(
    AFR = "#56B4E9", ADX = "#4B1FA8", EUR = "#fb8072"
  ),
  series.labels = c(
    Simulation_small = "Simulation small",
    Simulation_large = "Simulation large",
    Empirical = "Empirical"
  )
)


# internal functions ----


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
    add.population.roles()

  return(data)
}


# read chromosome-labelled kinship files for one source
read.kinship.chromosomes <- function(
    data.directory, chromosomes, data.type.input
) {
  paths <- file.path(
    path.expand(data.directory),
    paste0("kinship_unrelated.chr", chromosomes, ".parquet")
  )
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
  sizes <- candidates %>%
    count(data.type, rep, role, chrom, name = "available")
  if (any(sizes$available < downsample.size)) {
    stop(paste0(
      "A simulation kinship group contains fewer than ",
      downsample.size, " unique IDs"
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
    group_by(data.type, rep, role, chrom, sample.set) %>%
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
        fraction = histogram$counts / sum(histogram$counts)
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
      data.type, role, chrom, sample.set, xmin, xmax, xmid
    ) %>%
    summarise(
      mean.fraction = mean(fraction),
      sd.fraction = sd(fraction),
      replicate.count = n_distinct(rep),
      .groups = "drop"
    )
  empirical <- data %>%
    filter(data.type == "Empirical") %>%
    select(
      data.type, role, chrom, sample.set, xmin, xmax, xmid,
      mean.fraction = fraction
    ) %>%
    distinct() %>%
    mutate(sd.fraction = NA_real_, replicate.count = 1L)
  summary <- bind_rows(simulation, empirical) %>%
    mutate(
      role = factor(role, levels = c("AFR", "ADX", "EUR")),
      chrom = factor(
        chrom, levels = c(SELECTED.CHROMOSOMES, "all")
      ),
      data.type = factor(
        data.type,
        levels = c(
          "Simulation_small", "Simulation_large", "Empirical"
        )
      )
    )

  return(summary)
}


# construct the pairwise kinship distribution plot
make.kinship.plot <- function(data, plot.sample.set, breaks, styles) {
  plot.data <- data %>%
    filter(
      data.type == "Empirical" |
        (data.type != "Empirical" & sample.set == plot.sample.set)
    )
  dodge <- position_dodge(width = diff(breaks)[1] * 0.9)
  plot <- ggplot(
    plot.data,
    aes(
      x = xmid, y = mean.fraction, fill = role,
      group = role
    )
  ) +
    geom_col(
      position = dodge, width = diff(breaks)[1] * 0.85,
      color = "black", linewidth = 0.1
    ) +
    # geom_errorbar(
    #   data = plot.data,
    #   aes(
    #     ymin = pmax(0, mean.fraction - 2 * sd.fraction),
    #     ymax = mean.fraction + 2 * sd.fraction
    #   ),
    #   position = dodge, width = 0, linewidth = 0.45,
    #   color = "black", na.rm = TRUE
    # ) +
    facet_grid(
      chrom ~ data.type, drop = FALSE, scales = "free_y",
      labeller = labeller(
        data.type = as_labeller(styles$series.labels)
      )
    ) +
    scale_fill_manual(values = styles$population.colors) +
    labs(
      x = "Pairwise Kinship", y = "Fraction of pairs",
      title = "Pairwise Kinship Distributions Across Chromosomes",
      fill = NULL
    ) +
    xlim(-0.2, 0.0442) +
    guides(fill = guide_legend(order = 1)) +
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


# read selected-chromosome simulation and empirical kinship estimates
sim.small.kinship <- read.kinship.chromosomes(
  SIM.SMALL.DATA.DIR, SELECTED.CHROMOSOMES, "Simulation_small"
)
sim.large.kinship <- read.kinship.chromosomes(
  SIM.LARGE.DATA.DIR, SELECTED.CHROMOSOMES, "Simulation_large"
)
emp.chromosome.kinship <- read.kinship.chromosomes(
  EMPIRICAL.DATA.DIR, SELECTED.CHROMOSOMES, "Empirical"
)
emp.genome.kinship <- read.empirical.kinship.genome(
  EMPIRICAL.DATA.DIR
)

# build full and fixed-size downsampled simulation pair sets
simulation.kinship <- bind_rows(
  sim.small.kinship, sim.large.kinship
)
selected.kinship.ids <- select.kinship.ids(
  simulation.kinship, DOWNSAMPLE.SIZE, RANDOM.SEED
)
simulation.kinship <- apply.kinship.selection(
  simulation.kinship, selected.kinship.ids
)
empirical.kinship <- bind_rows(
  emp.chromosome.kinship, emp.genome.kinship
) %>%
  mutate(sample.set = "full")

# summarize common-bin histograms and construct the primary plot
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
kinship.plot <- make.kinship.plot(
  kinship.summary, PLOT.SAMPLE.SET, kinship.breaks, PLOT.STYLES
)
print(kinship.plot)

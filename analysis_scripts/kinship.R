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


SIM.SMALL.DATA.DIR <- "~/scratch/OOA_NAAdmixture_small/stats"
SIMDOWN.SMALL.DATA.DIR <- "~/scratch/OOA_NAAdmixture_smallOnekgDownsample/stats"
SIM.LARGE.DATA.DIR <- "~/scratch/OOA_NAAdmixture_large/stats"
SIMDOWN.LARGE.DATA.DIR <- "~/scratch/OOA_NAAdmixture_largeOnekgDownsample/stats"
EMPIRICAL.DATA.DIR <- "~/scratch/OOA_NAAdmixture_1kG/stats"
CHROMOSOMES <- as.character(1:22)
SELECTED.CHROMOSOMES <- c("1")
SOURCE.LEVELS <- c(
  "Simulation_small", "Simulation_large", "Simulation_small_simDown",
  "Simulation_large_simDown", "Empirical"
)
KINSHIP.BIN.WIDTH <- 0.01
PLOT.BASE.SIZE <- 24
PLOT.STYLES <- list(
  population.colors = c(
    AFR = "#56B4E9", ADX = "#4B1FA8", EUR = "#fb8072",
    YRI = "#eec4dc", ASW = "#e44b8d", CEU = "#bb437e"
  ),
  source.colors = c(
    Simulation_small = "#9A83CE", Simulation_large = "#32146F",
    Simulation_small_simDown = "#6F55B5",
    Simulation_large_simDown = "#4B1FA8"
  ),
  series.labels = c(
    Simulation_small = "Sm. Sim.",
    Simulation_small_simDown = "Sm. D. Sim.",
    Simulation_large = "Lg. Sim.",
    Simulation_large_simDown = "Lg. D. Sim.",
    Empirical = "Emp."
  )
)


# internal functions ----


# describe the kinship estimator, scope, and uncertainty concisely
kinship.plot.subtitle <- function(chromosomes) {
  subtitle <- paste0(
    "KING estimator · chromosomes ", paste(chromosomes, collapse = ", "),
    " and empirical genome-wide · ±2 SD simulations"
  )
  return(subtitle)
}


# retain source-specific populations before any histogram calculations
apply.kinship.source.contract <- function(data) {
  retained <- data %>%
    filter(
      (data.type %in% c(
        "Simulation_small", "Simulation_small_simDown"
      ) & pop %in% c("AFR", "ADX", "EUR")) |
        (data.type %in% c(
          "Simulation_large", "Simulation_large_simDown"
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
    glue("kinship_unrelated.chr{chromosomes}.parquet")
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


# construct the pairwise kinship distribution plot
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
    stop("Unsupported kinship plot view: ", view)
  )
  filtered <- data %>%
    filter(as.character(data.type) %in% sources) %>%
    filter(view != "simulated" | pop == "ADX") %>%
    mutate(data.type = factor(as.character(data.type), levels = sources)) %>%
    arrange(data.type)
  return(filtered)
}


# construct one scoped pairwise kinship distribution plot
make.kinship.plot <- function(data, breaks, styles, view) {
  plot.data <- filter.plot.view(data, view) %>%
    mutate(chrom = forcats::fct_drop(as.factor(chrom))) %>%
    mutate(plot.key = if (view == "simulated") {
      as.character(data.type)
    } else {
      as.character(pop)
    })
  chromosome.scope <- plot.data %>%
    filter(data.type != "Empirical") %>%
    pull(chrom) %>%
    as.character() %>%
    unique()
  subtitle <- kinship.plot.subtitle(chromosome.scope)
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
    scale_fill_manual(values = if (view == "simulated") {
      styles$source.colors
    } else {
      styles$population.colors
    }, labels = if (view == "simulated") styles$series.labels else waiver()) +
    labs(
      x = "Pairwise Kinship", y = "Fraction of pairs",
      title = if (view == "simulated") {
        "Pairwise Kinship Distributions: Simulated ADX"
      } else {
        "Pairwise Kinship Distributions: Small Simulation and Empirical"
      },
      subtitle = subtitle,
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
simDown.small.kinship <- read.kinship.chromosomes(
  SIMDOWN.SMALL.DATA.DIR, SELECTED.CHROMOSOMES,
  "Simulation_small_simDown"
)
sim.large.kinship <- read.kinship.chromosomes(
  SIM.LARGE.DATA.DIR, SELECTED.CHROMOSOMES, "Simulation_large"
)
simDown.large.kinship <- read.kinship.chromosomes(
  SIMDOWN.LARGE.DATA.DIR, SELECTED.CHROMOSOMES,
  "Simulation_large_simDown"
)
emp.chromosome.kinship <- read.kinship.chromosomes(
  EMPIRICAL.DATA.DIR, SELECTED.CHROMOSOMES, "Empirical"
)
emp.genome.kinship <- read.empirical.kinship.genome(
  EMPIRICAL.DATA.DIR
)

# combine all unrelated simulation and empirical pairs
simulation.kinship <- bind_rows(
  sim.small.kinship, simDown.small.kinship,
  sim.large.kinship, simDown.large.kinship
)
empirical.kinship <- bind_rows(
  emp.chromosome.kinship, emp.genome.kinship
)

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
kinship.small.empirical.plot <- make.kinship.plot(
  kinship.summary, kinship.breaks, PLOT.STYLES, "small_empirical"
)
kinship.simulated.plot <- make.kinship.plot(
  kinship.summary, kinship.breaks, PLOT.STYLES, "simulated"
)
print(kinship.small.empirical.plot)
print(kinship.simulated.plot)

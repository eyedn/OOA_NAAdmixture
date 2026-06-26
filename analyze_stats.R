# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology 
# Mooney Lab
# ___
# analyze_stats.R
# ______________________________________________________________________________


library(tidyverse)
library(nanoparquet)
library(pillar)
library(glue)


# data paths ----
data.dir <- "~/scratch/OOA_testing/OOA_NAAdmixture_chr1/stats"
ancestry.file <- file.path(data.dir, "ancestry.parquet")
pi.theta.file <- file.path(data.dir, "pi_theta_stats.parquet")
sfs.file <- file.path(data.dir, "sfs.parquet")
sfs.2d.file <- file.path(data.dir, "sfs_2d.parquet")
ld.decay.file <- file.path(data.dir, "ld_decay.parquet")
kinship.file <- file.path(data.dir, "kinship.parquet")

# ancestry ----
read_parquet_info(ancestry.file)
read_parquet_schema(ancestry.file)
ancestry <- read_parquet(ancestry.file)

summary(ancestry %>% filter(pop == "ADX") %>% pull(afr_tspop))

admix.plot <- ancestry %>%
  mutate(pop = factor(pop, levels = c("AFR", "ADX", "EUR"))) %>%
  arrange(pop, desc(afr_tspop)) %>%
  mutate(ind = row_number())

boundaries <- admix.plot %>%
  count(pop) %>%
  mutate(
    xmax = cumsum(n),
    xmin = lag(xmax, default = 0) + 1,
    midpoint = (xmin + xmax) / 2
  )

admix.long <- admix.plot %>%
  pivot_longer(
    c(afr_q, eur_q),
    names_to = "ancestry",
    values_to = "q"
  )

ggplot(admix.long,
       aes(ind, q, fill = ancestry)) +
  geom_col(width = 1) +
  geom_vline(
    xintercept = boundaries$xmax[-nrow(boundaries)] + 0.5,
    linewidth = 0.4
  ) +
  scale_fill_manual(
    values = c(
      afr_q = "#56B4E9",
      eur_q = "#fb8072"
    )
  ) +
  scale_x_continuous(
    breaks = boundaries$midpoint,
    labels = boundaries$pop
  ) +
  labs(x = NULL, y = "Ancestry proportion", fill = NULL) +
  theme_bw() +
  theme(
    axis.ticks.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )

# diversity
read_parquet_info(pi.theta.file)
read_parquet_schema(pi.theta.file)
pi.theta <- read_parquet(pi.theta.file)

# sfs ----
read_parquet_info(sfs.file)
read_parquet_schema(sfs.file)
sfs <- read_parquet(sfs.file)

pop_cols <- c(
  AFR = "#56B4E9",
  ADX = "#4B1FA8",
  EUR = "#fb8072"
)

sfs.plot <- sfs %>%
  filter(
    derived_allele_count > 0,
    derived_allele_count <= 15
  ) %>%
  mutate(
    pop = factor(pop, levels = c("AFR", "ADX", "EUR"))
  ) %>%
  group_by(pop) %>%
  mutate(prop = count / sum(count)) %>%
  ungroup()

pd <- position_dodge(width = 0.8)

ggplot(
  sfs.plot,
  aes(
    x = factor(derived_allele_count),
    y = prop,
    fill = pop,
    color = pop,
    group = pop
  )
) +
  geom_col(
    position = pd,
    width = 0.75
  ) +
  scale_fill_manual(values = pop_cols) +
  scale_color_manual(values = pop_cols) +
  labs(
    x = "Derived allele count",
    y = "Proportion of segregating sites",
    fill = "Population",
    color = "Population"
  ) +
  theme_bw() +
  theme(
    legend.position = "top"
  )

# 2d sfs ----
read_parquet_info(sfs.2d.file)
read_parquet_schema(sfs.2d.file)
sfs.2d <- read_parquet(sfs.2d.file)

sfs2d.plot <- sfs.2d %>%
  filter(
    pop1_count <= 50,
    pop2_count <= 50
  ) %>%
  mutate(
    pair = factor(
      paste(pop1, pop2, sep = " vs "),
      levels = c(
        "AFR vs EUR",
        "AFR vs ADX",
        "EUR vs ADX"
      )
    )
  )

ggplot(
  sfs2d.plot,
  aes(
    pop1_count,
    pop2_count,
    fill = count
  )
) +
  geom_tile() +
  scale_fill_gradientn(
    colours = c("white", "gold", "orange", "red", "darkred"),
    trans = "log10"
  ) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
  facet_wrap(~pair, nrow = 1) +
  scale_x_continuous(breaks = seq(0, 50, 5)) +
  scale_y_continuous(breaks = seq(0, 50, 5)) +
  coord_equal() +
  labs(
    x = "Derived allele count in population 1",
    y = "Derived allele count in population 2"
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.position = "right"
  )

# ld ----
read_parquet_info(ld.decay.file)
read_parquet_schema(ld.decay.file)
ld.decay <- read_parquet(ld.decay.file)

# look at one replicate
rep.use <- 50
min.bin <- 20000
max.bin <- 200000

ld.plot <- ld.decay %>%
  filter(
    rep == rep.use,
    distance_bin_bp >= min.bin,
    distance_bin_bp <= max.bin
  ) %>%
  group_by(pop, distance_bin_bp) %>%
  summarise(
    sum_r2 = sum(sum_r2),
    n_pairs = sum(n_pairs),
    mean_r2 = sum_r2 / n_pairs,
    .groups = "drop"
  ) %>%
  mutate(
    pop = factor(pop, levels = c("AFR", "ADX", "EUR"))
  )

ggplot(
  ld.plot,
  aes(
    x = distance_bin_bp,
    y = mean_r2,
    color = pop
  )
) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = pop_cols) +
  labs(
    x = "Distance between SNPs (bp)",
    y = expression(r^2),
    color = "Population"
  ) +
  theme_bw() +
  theme(
    legend.position = "top"
  )

# now look at all replicates

min.bin <- 20000
max.bin <- 200000

# aggregate weighted mean of r2 across windows within each replicate
ld.rep <- ld.decay %>%
  filter(
    distance_bin_bp >= min.bin,
    distance_bin_bp <= max.bin
  ) %>%
  group_by(rep, pop, distance_bin_bp) %>%
  summarise(
    sum_r2 = sum(sum_r2),
    n_pairs = sum(n_pairs),
    av_r2 = sum_r2 / n_pairs,
    .groups = "drop"
  )

# mean and SE of mean r2 across replicates
ld.plot <- ld.rep %>%
  group_by(pop, distance_bin_bp) %>%
  summarise(
    mean_r2 = mean(av_r2),
    sd_r2   = sd(av_r2),
    se_r2   = sd_r2 / sqrt(n()),
    lower   = mean_r2 - 1.96 * se_r2,
    upper   = mean_r2 + 1.96 * se_r2,
    .groups = "drop"
  ) %>%
  mutate(
    pop = factor(pop, levels = c("AFR", "ADX", "EUR"))
  )

ggplot(
  ld.plot,
  aes(
    x = distance_bin_bp,
    y = mean_r2,
    color = pop,
    fill = pop
  )
) +
  geom_ribbon(
    aes(ymin = lower, ymax = upper),
    alpha = 0.2,
    color = NA
  ) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = pop_cols) +
  scale_fill_manual(values = pop_cols) +
  labs(
    x = "Distance between SNPs (bp)",
    y = expression(r^2),
    color = "Population",
    fill = "Population"
  ) +
  theme_bw() +
  theme(
    legend.position = "top"
  )


# kinship ----
read_parquet_info(kinship.file)
read_parquet_schema(kinship.file)
kinship <- read_parquet(kinship.file)

kinship %>%
  filter(rep == 5) %>%
  ggplot(aes(x = kinship, fill = pop, color = pop)) +
  geom_density(alpha = 0.25, linewidth = 1) +
  scale_fill_manual(values = pop_cols) +
  scale_color_manual(values = pop_cols) +
  theme_bw()

kinship %>%
  filter(rep == 1) %>%
  ggplot(aes(pop, kinship, fill = pop, color = pop)) +
  geom_violin(alpha = 0.3) +
  geom_boxplot(width = 0.15, fill = "white") +
  scale_fill_manual(values = pop_cols) +
  scale_color_manual(values = pop_cols) +
  theme_bw()

kinship %>%
  filter(rep == 1) %>%
  ggplot(aes(kinship, color = pop)) +
  stat_ecdf(linewidth = 1) +
  scale_color_manual(values = pop_cols) +
  theme_bw()

# idea 1
# Common x-grid for all density estimates
x_grid <- seq(
  min(kinship$kinship),
  max(kinship$kinship),
  length.out = 500
)

# Estimate density separately for each population and replicate
densities <- kinship %>%
  group_by(pop, rep) %>%
  summarise(
    density = list(
      density(
        kinship,
        from = min(x_grid),
        to = max(x_grid),
        n = length(x_grid)
      )
    ),
    .groups = "drop"
  ) %>%
  mutate(
    x = map(density, "x"),
    y = map(density, "y")
  ) %>%
  select(-density) %>%
  unnest(c(x, y))

# Mean ± SD across replicates
density.summary <- densities %>%
  group_by(pop, x) %>%
  summarise(
    mean_density = mean(y),
    sd_density = sd(y),
    lower = pmax(mean_density - sd_density, 0),
    upper = mean_density + sd_density,
    .groups = "drop"
  ) %>%
  mutate(
    pop = factor(pop, levels = c("AFR", "ADX", "EUR"))
  )

ggplot(
  density.summary,
  aes(
    x = x,
    y = mean_density,
    color = pop,
    fill = pop
  )
) +
  geom_ribbon(
    aes(ymin = lower, ymax = upper),
    alpha = 0.20,
    color = NA
  ) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = pop_cols) +
  scale_fill_manual(values = pop_cols) +
  labs(
    x = "Kinship coefficient",
    y = "Density",
    color = "Population",
    fill = "Population"
  ) +
  theme_bw() +
  theme(
    legend.position = "top"
  )

# idea 2
kinship.summary <- kinship %>%
  group_by(rep, pop) %>%
  summarise(
    mean_kinship = mean(kinship),
    median_kinship = median(kinship),
    sd_kinship = sd(kinship),
    .groups = "drop"
  ) %>%
  mutate(
    pop = factor(pop, levels = c("AFR", "ADX", "EUR"))
  )
ggplot(
  kinship.summary,
  aes(
    x = pop,
    y = mean_kinship,
    fill = pop,
    color = pop
  )
) +
  geom_violin(alpha = 0.25, linewidth = 1) +
  geom_boxplot(
    width = 0.15,
    fill = "white",
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.08,
    size = 2,
    alpha = 0.6
  ) +
  scale_fill_manual(values = pop_cols) +
  scale_color_manual(values = pop_cols) +
  labs(
    x = "Population",
    y = "Mean pairwise kinship"
  ) +
  theme_bw() +
  theme(
    legend.position = "none"
  )

kinship.summary2 <- kinship.summary %>%
  group_by(pop) %>%
  summarise(
    mean = mean(mean_kinship),
    sd = sd(mean_kinship),
    .groups = "drop"
  )

ggplot(
  kinship.summary2,
  aes(
    pop,
    mean,
    fill = pop,
    color = pop
  )
) +
  geom_col(width = 0.6, alpha = 0.8) +
  geom_errorbar(
    aes(
      ymin = mean - sd,
      ymax = mean + sd
    ),
    width = 0.15,
    linewidth = 1
  ) +
  scale_fill_manual(values = pop_cols) +
  scale_color_manual(values = pop_cols) +
  theme_bw() +
  theme(
    legend.position = "none"
  )

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


# read data ----
data.dir <- "~/scratch/OOA_NAAdmixture_chr22/stats"
admix.q.file <- file.path(data.dir, "admixture_q.parquet")
pi.theta.file <- file.path(data.dir, "pi_theta_stats.parquet")
sfs.file <- file.path(data.dir, "sfs.parquet")
sfs.2d.file <- file.path(data.dir, "sfs_2d.parquet")
ld.decay.file <- file.path(data.dir, "ld_decay.parquet")
kinship.file <- file.path(data.dir, "kinship.parquet")

# admixture Q matrix inference -----
read_parquet_info(admix.q.file)
read_parquet_schema(admix.q.file)
admix.q <- read_parquet(admix.q.file)

# pi and theta -----
read_parquet_info(pi.theta.file)
read_parquet_schema(pi.theta.file)
pi.theta <- read_parquet(pi.theta.file)

# 1D sfs -----
read_parquet_info(sfs.file)
read_parquet_schema(sfs.file)
sfs <- read_parquet(sfs.file)

# 2D sfs -----
read_parquet_info(sfs.2d.file)
read_parquet_schema(sfs.2d.file)
sfs.2d <- read_parquet(sfs.2d.file)

# ld decay -----
read_parquet_info(ld.decay.file)
read_parquet_schema(ld.decay.file)
ld.decay <- read_parquet(ld.decay.file)

# KING kinship coefficients -----
read_parquet_info(kinship.file)
read_parquet_schema(kinship.file)
kinship <- read_parquet(kinship.file)

# summarize admixture Q matrix results ----
summary(admix.q %>% filter(pop == "ADX") %>% pull(afr_q))

admix.plot <- admix.q %>%
  mutate(pop = factor(pop, levels = c("AFR", "ADX", "EUR"))) %>%
  arrange(pop, desc(afr_q)) %>%
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

# sfs ----
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
    width = 0.75,
    alpha = 0.8
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


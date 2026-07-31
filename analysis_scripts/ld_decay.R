# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology 
# Mooney Lab
# ___
# ld_decay.R
# ______________________________________________________________________________


library(tidyverse)
library(nanoparquet)
library(glue)
library(ggridges)


sim.data.dir <- "~/scratch/OOA_NAAdmixture_small_cp/stats"
emp.data.dir <- "~/scratch/OOA_NAAdmixture_1kG/stats"


ld.xlim.lower <- 0
ld.xlim.upper <- 250000
ld.xbreaks.freq <- 50000
ld.xbreaks <- seq(ld.xlim.lower, ld.xlim.upper, ld.xbreaks.freq)

# read in data ----


# read in simulated data ----
sim.ld.chr <- map_dfr(
  1:22,
  \(chr) {
    file.path(sim.data.dir, glue("ld_decay.chr{chr}.parquet")) |>
      read_parquet() %>%
      mutate(
        chrom = "all", rep = as.numeric(rep), data.type = "Simulation",
        rep = as.numeric(rep), distance_bin_bp = as.numeric(distance_bin_bp),
        mean_r2 = as.numeric(mean_r2), sum_r2 = as.numeric(sum_r2),
        n_pairs = as.numeric(n_pairs) 
      )
    }
  )

sim.ld.genome <- read_parquet(
  file.path(sim.data.dir, "ld_decay.parquet")
  ) %>%
  mutate(
    chrom = "all", data.type = "Simulation", rep = as.numeric(rep),
    distance_bin_bp = as.numeric(distance_bin_bp), 
    mean_r2 = as.numeric(mean_r2), sum_r2 = as.numeric(sum_r2),
    n_pairs = as.numeric(n_pairs)
    )


# read in empirical data ----
emp.ld.chr <- map_dfr(
  1:22,
  \(chr) {
    file.path(emp.data.dir, glue("ld_decay.chr{chr}.parquet")) |>
      read_parquet() %>%
      mutate(chrom = "all", data.type = "Empirical")
    }
  )

emp.ld.genome <- read_parquet(
  file.path(emp.data.dir, "ld_decay.parquet")
  ) %>%
  mutate(chrom = "all", data.type = "Empirical")


## combine data ----
ld.genome <- rbind(sim.ld.genome, emp.ld.genome) %>%
  mutate(
    role = ifelse(
      pop %in% c("AFR", "YRI"), "AFR",
      ifelse(pop %in% c("EUR", "CEU"), "EUR", "ADX")
    )
  )

# summarize ld decay ----
ld.summary <- ld.genome %>%
  mutate(
    pop = factor(pop, levels = c("AFR", "ADX", "EUR", "YRI", "ASW", "CEU"))
  ) %>%
  group_by(pop, role, data.type, distance_bin_bp) %>%
  summarize(
    mean = mean(mean_r2, na.rm = TRUE),
    sd   = sd(mean_r2, na.rm = TRUE),
    .groups = "drop"
  )


## plot ld decay, facet by data.type ----
ld.facetDataType.plot <- ggplot(
  ld.summary,
  aes(x = distance_bin_bp, y = mean, color = pop, fill = pop)
  ) +
  geom_ribbon(
    data = ld.summary %>% filter(data.type == "Simulation"),
    aes(ymin = pmax(0, mean - 2 * sd), ymax = mean + 2 * sd),
    alpha = 0.20, color = NA
    ) +
  geom_line(linewidth = 1.25) +
  scale_fill_manual(
    name = NULL,
    breaks = c("AFR", "ADX", "EUR", "YRI",  "ASW", "CEU"),
    values = c(
      "AFR" = "#56B4E9", "ADX" = "#4B1FA8", "EUR" = "#fb8072",
      "YRI" = "#eec4dc", "ASW" = "#e44b8d", "CEU" = "#bb437e"
      )
    ) +
  scale_color_manual(
    name = NULL,
    breaks = c("AFR", "ADX", "EUR", "YRI",  "ASW", "CEU"),
    values = c(
      "AFR" = "#56B4E9", "ADX" = "#4B1FA8", "EUR" = "#fb8072",
      "YRI" = "#eec4dc", "ASW" = "#e44b8d", "CEU" = "#bb437e"
      )
    ) +
  scale_x_continuous(
    limits = c(ld.xlim.lower, ld.xlim.upper), breaks = ld.xbreaks
    ) +
  facet_wrap(~data.type) +
  labs(
    x = "Distance between SNPs (bp)", y = expression("Mean " * r^2),
    title = "Genome-wide LD Decay across Replicates",
    subtitle = "Lines represent the mean LD decay across replicates; ribbons denote ±2 SD"
  ) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
  theme_bw(base_size = 24) +
  theme(
    legend.position = "top", legend.title = element_blank(), 
    panel.grid.minor = element_blank(), 
    strip.background = element_blank(), strip.text = element_text(face = "bold")
  )
ld.facetDataType.plot


## plot ld decay, facet by role ----
ld.facetRole.plot <- ggplot(
  ld.summary,
  aes(x = distance_bin_bp, y = mean, color = pop, fill = pop)
  ) +
  geom_ribbon(
    data = ld.summary %>% filter(data.type == "Simulation"),
    aes(ymin = pmax(0, mean - 2 * sd), ymax = mean + 2 * sd),
    alpha = 0.20, color = NA
  ) +
  geom_line(linewidth = 1.25) +
  scale_fill_manual(
    name = NULL,
    breaks = c("AFR", "ADX", "EUR", "YRI",  "ASW", "CEU"),
    values = c(
      "AFR" = "#56B4E9", "ADX" = "#4B1FA8", "EUR" = "#fb8072",
      "YRI" = "#eec4dc", "ASW" = "#e44b8d", "CEU" = "#bb437e"
    )
  ) +
  scale_color_manual(
    name = NULL,
    breaks = c("AFR", "ADX", "EUR", "YRI",  "ASW", "CEU"),
    values = c(
      "AFR" = "#56B4E9", "ADX" = "#4B1FA8", "EUR" = "#fb8072",
      "YRI" = "#eec4dc", "ASW" = "#e44b8d", "CEU" = "#bb437e"
    )
  ) +
  scale_x_continuous(
    limits = c(ld.xlim.lower, ld.xlim.upper), breaks = ld.xbreaks
  ) +
  facet_wrap(~role) +
  labs(
    x = "Distance between SNPs (bp)", y = expression("Mean " * r^2),
    title = "Genome-wide LD Decay across Replicates",
    subtitle = "Lines represent the mean LD decay across replicates; ribbons denote ±2 SD"
  ) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
  theme_bw(base_size = 24) +
  theme(
    legend.position = "top", legend.title = element_blank(), 
    panel.grid.minor = element_blank(),
    strip.background = element_blank(), strip.text = element_text(face = "bold")
  )
ld.facetRole.plot


## plot ld decay, all population/data.types ----
ld.noRibbon.plot <- ggplot(
  ld.summary,
  aes(x = distance_bin_bp, y = mean, color = pop)
) +
  geom_line(
    linewidth = 2,
    aes(linetype = data.type)
  ) +
  scale_color_manual(
    name = NULL,
    breaks = c("AFR", "ADX", "EUR", "YRI", "ASW", "CEU"),
    values = c(
      "AFR" = "#56B4E9",
      "ADX" = "#4B1FA8",
      "EUR" = "#fb8072",
      "YRI" = "#eec4dc",
      "ASW" = "#e44b8d",
      "CEU" = "#bb437e"
    )
  ) +
  scale_linetype_manual(
    name = NULL,
    breaks = c("Simulation", "Empirical"),
    values = c(
      "Simulation" = "dashed",
      "Empirical" = "solid"
    )
  ) +
  scale_x_continuous(
    limits = c(ld.xlim.lower, ld.xlim.upper),
    breaks = ld.xbreaks
  ) +
  labs(
    x = "Distance between SNPs (bp)",
    y = expression("Mean " * r^2),
    title = "Genome-wide LD Decay across Replicates"
  ) +
  guides(
    color = guide_legend(
      order = 1,
      nrow = 1,
      byrow = TRUE,
      override.aes = list(linetype = "solid")
    ),
    linetype = guide_legend(
      order = 2,
      override.aes = list(
        color = "black",
        linewidth = 1.25
      )
    )
  ) +
  theme_bw(base_size = 24) +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    legend.key.width = unit(42, "pt"),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  )
ld.noRibbon.plot

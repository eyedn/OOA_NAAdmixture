# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology 
# Mooney Lab
# ___
# diversity.R
# ______________________________________________________________________________


# set up ----
library(tidyverse)
library(nanoparquet)
library(glue)
library(ggridges)
library(broom)
library(car)
library(emmeans)


sim.data.dir <- "~/scratch/OOA_NAAdmixture_small/stats"
emp.data.dir <- "~/scratch/OOA_NAAdmixture_1kG/stats"

pd.width <- 0.75
pd <- position_dodge(width = pd.width)


# read in data ----


## read in simulated data ----
sim.div.chr <- map_dfr(
  1:22,
  \(chr) {
    file.path(sim.data.dir, glue("pi_theta_stats.chr{chr}.parquet")) |>
      read_parquet() |>
      mutate(
        rep = as.numeric(rep),mutation_rate = as.numeric(mutation_rate),
        wattersons_const = as.numeric(wattersons_const),
        segregating_sites = as.numeric(segregating_sites),
        ne_value = as.numeric(ne_value), value = as.numeric(value), 
        data.type = "sim"
        )
    }
  )

sim.div.genome <- read_parquet(
  file.path(sim.data.dir, glue("pi_theta_stats.parquet"))
  ) %>%
  mutate(
    rep = as.numeric(rep), mutation_rate = as.numeric(mutation_rate),
    span = NA, wattersons_const = as.numeric(wattersons_const),
    segregating_sites = as.numeric(segregating_sites),
    ne_value = as.numeric(ne_value), chrom = "all", value = as.numeric(value), 
    data.type = "sim"
    )

## read in empirical intergenic data ----
emp.div.intergenic.chr <- map_dfr(
  1:22,
  \(chr) {
    file.path(
      emp.data.dir, glue("pi_theta_stats_intergenic.chr{chr}.parquet")
      ) %>%
      read_parquet() %>%
      mutate(
        chrom = as.character(chrom), span = as.numeric(span),
        value = as.numeric(value), data.type = "intergenic"
      )
    }
  )

emp.div.intergenic.genome <- read_parquet(
  file.path(emp.data.dir, glue("pi_theta_stats_intergenic.parquet"))
  ) %>%
  mutate(
    chrom = "all", span = NA, value = as.numeric(value), 
    data.type = "intergenic"
    )

## read in empirical full-callable data ----
emp.div.fullCall.chr <- map_dfr(
  1:22,
  \(chr) {
    file.path(
      emp.data.dir, glue("pi_theta_stats_full_callable_chrom.chr{chr}.parquet")
    ) %>%
      read_parquet() %>%
      mutate(
        chrom = as.character(chrom), span = as.numeric(span),
        value = as.numeric(value),  data.type = "fullCall"
      )
    }
  )

emp.div.fullCall.genome <- read_parquet(
  file.path(emp.data.dir, glue("pi_theta_stats_full_callable_chrom.parquet"))
  ) %>%
  mutate(
    chrom = "all", span = NA, value = as.numeric(value), data.type = "fullCall"
    )


## combine diversity data ----
div <- rbind(
  sim.div.chr, sim.div.genome, emp.div.intergenic.chr,
  emp.div.intergenic.genome, emp.div.fullCall.chr, emp.div.fullCall.genome
) %>%
  mutate(
    role = ifelse(
      pop %in% c("AFR", "YRI"), "AFR",
      ifelse(pop %in% c("EUR", "CEU"), "EUR", "ADX")
    )
  )


# summarize data ----
div.sim.summary <- div %>%
  filter(stat %in% c("pi", "theta"), data.type == "sim") %>%
  group_by(stat, chrom, pop) %>%
  summarize(mean = mean(value), sd = sd(value), .groups = "drop") %>%
  mutate(pop = factor(pop, levels = c("AFR", "ADX", "EUR"))) %>%
  mutate(
    stat = factor(
      stat, levels = c("pi", "theta"),
      labels = c(expression(pi),expression(theta[w]))
      )
    )

div.emp.summary <- div %>%
  filter(stat %in% c("pi", "theta"), data.type != "sim") %>%
  mutate(
    dataset = factor(
      if_else(data.type == "intergenic", "Intergenic", "Full callable"),
      levels = c("Intergenic", "Full callable")
      ),
    pop = factor(
      pop, levels = c("YRI", "ASW", "CEU")
      ),
    stat = factor(
      stat, levels = c("pi", "theta"), 
      labels = c(expression(pi), expression(theta[w]))
      )
    )

div.analysis <- bind_rows(
  div.sim.summary %>%
    mutate(span = NA) %>%
    transmute(
      stat, chrom, span, dataset = "Simulation", population = pop,
      estimate = mean
      ),
  div.emp.summary %>%
    transmute(
      stat, chrom, span, dataset = as.character(dataset), population = pop,
      estimate = value
      )
  ) %>%
  left_join(
    chrom.lengths %>% transmute(chrom, chr.len.mb = chr_len / 1e6), by = "chrom"
    ) %>%
  mutate(
    chr.len.mb = if_else(dataset == "Simulation", chr.len.mb, as.numeric(span))
    ) %>%
  filter(chrom != "all") %>%
  select(-span)


## stats ----


### ancova ----
div.models <- div.analysis %>%
  group_by(stat, dataset) %>%
  nest()

div.models <- div.models %>%
  mutate(fit = map(data, ~ lm(estimate ~ population + chr.len.mb,data = .x)))


### anova ----
div.models <- div.models %>% mutate(anova = map(fit, ~ Anova(.x, type = 2)))

div.anova.table <- div.models %>%
  transmute(stat, dataset, anova = map(anova, tidy)) %>%
  unnest(anova) %>%
  filter(term != "Residuals")
div.anova.table


### compute estiamted marginal means ----
div.models <- div.models %>%
  mutate(emmeans = map(fit, ~ emmeans(.x, ~ population)))

div.emmeans.table <- div.models %>%
  transmute(stat, dataset, emmeans = map(emmeans, tidy)) %>%
  unnest(emmeans)
div.emmeans.table

div.models <- div.models %>%
  mutate(
    tukey = map(
      emmeans, ~ summary(
        pairs(.x, adjust = "tukey"), infer = c(TRUE, TRUE)
        ) %>% 
        as_tibble()
      )
    )

div.tukey.table <- div.models %>%
  select(stat, dataset, tukey) %>%
  unnest(tukey)
div.tukey.table


### compute difference and fold change table ----
div.summary.table <- div.tukey.table %>%
  separate(contrast, into = c("pop1", "pop2"), sep = " - ")

div.summary.table <- div.summary.table %>%
  left_join(
    div.emmeans.table %>%
      select(stat, dataset, pop1 = population, mean1 = estimate),
    by = c("stat", "dataset", "pop1")
    ) %>%
  left_join(
    div.emmeans.table %>%
      select(stat, dataset, pop2 = population, mean2 = estimate),
    by = c("stat", "dataset", "pop2")
    ) %>%
  mutate(
    fold.change = mean1 / mean2,
    percent.change = 100 * (mean1 - mean2) / mean2
  ) %>%
  transmute(
    stat, dataset, comparison = paste(pop1, "vs", pop2), 
    adjusted.mean.1 = mean1, adjusted.mean.2 = mean2, difference = estimate,
    fold.change, percent.change, lower.CL, upper.CL, p.value
  )
div.summary.table


## dot plot ----
chrom.use <- c("all", as.character(1:22))
div.plot <- ggplot(
  div.sim.summary %>% filter(chrom %in% chrom.use),
  aes(x = factor(chrom, levels = chrom.levels), y = mean)
  ) +
  geom_errorbar(
    aes(ymin = mean - 2 * sd, ymax = mean + 2 * sd, group = pop), 
    color = "black", position = pd, width = 0.15, linewidth = 0.9
    ) +
  geom_point(
    aes(fill = pop, shape = "Simulation", group = pop), position = pd,
    color = "black", stroke = 1, size = 3.5
    ) +
  geom_point(
    data = div.emp.summary %>% filter(chrom %in% chrom.use),
    aes(
      x = factor(chrom, levels = chrom.levels), y = value, fill = pop, 
      shape = dataset, group = pop
      ),
    inherit.aes = FALSE, position = pd, color = "black", stroke = 1, size = 3.5
    ) +
  facet_wrap(~stat, scales = "free_y", labeller = label_parsed, ncol = 1) +
  scale_fill_manual(
    name = NULL,
    breaks = c("AFR", "YRI", "ADX", "ASW", "EUR", "CEU"),
    values = c(
      "AFR" = "#56B4E9", "ADX" = "#4B1FA8", "EUR" = "#fb8072",
      "YRI" = "#eec4dc", "ASW" = "#e44b8d", "CEU" = "#bb437e"
    )
  ) +
  scale_shape_manual(
    name = NULL,
    breaks = c("Simulation", "Intergenic", "Full callable"),
    values = c("Simulation" = 21, "Intergenic" = 24, "Full callable" = 23)
    ) +
  labs(
    x = "Chromosome", y = NULL, title = "Genetic diversity statistics",
    subtitle = "Error bars represent ±2 SD"
    ) +
  guides(
    fill = guide_legend(
      order = 1, override.aes = list(shape = 21, color = "black", size = 3.5)
      ),
    shape = guide_legend(
      order = 2, 
      override.aes = list(fill = "white", color = "black", size = 3.5)
      )
    ) +
  theme_bw(base_size = 24) +
  theme(
    legend.position = "top", legend.direction = "horizontal",
    legend.title = element_blank(), panel.grid.minor = element_blank(),
    strip.background = element_blank(), 
    strip.text = element_text(face = "bold"), panel.spacing = unit(1, "lines")
  )
div.plot


## chrom hist plot ----
chrom.div.box.df <- bind_rows(
  div.sim.summary %>%
    filter(chrom != "all") %>%
    transmute(stat, chrom, dataset = "Simulation", pop, estimate = mean),
  div.emp.summary %>%
    filter(chrom != "all") %>%
    transmute(stat, chrom, dataset, pop, estimate = value)
  ) %>%
  mutate(
    dataset = factor(
      dataset, levels = c("Simulation", "Intergenic", "Full callable")
      ),
    pop = factor(
      pop, levels = c("AFR", "YRI", "ADX", "ASW", "EUR", "CEU")
      )
    )

chrom.div.box.plot <- ggplot(
  chrom.div.box.df, aes(x = dataset, y = estimate, fill = pop)
  ) +
  geom_boxplot(
    aes(group = interaction(dataset, pop)), position = pd, width = 0.6,
    outlier.shape = NA, color = "black", linewidth = 0.8
    ) +
  geom_jitter(
    aes(group = interaction(dataset, pop)),
    position = position_jitterdodge(
      dodge.width = pd.width, jitter.width = 0.08
      ),
    shape = 21, color = "black", stroke = 0.7, size = 2.5
    ) +
  facet_wrap(~stat, scales = "free_y", labeller = label_parsed, ncol = 2) +
  scale_fill_manual(
    name = NULL, breaks = c("AFR", "ADX", "EUR", "YRI", "ASW", "CEU"),
    values = c(
      "AFR" = "#56B4E9", "ADX" = "#4B1FA8", "EUR" = "#fb8072",
      "YRI" = "#eec4dc", "ASW" = "#e44b8d", "CEU" = "#bb437e"
      )
    ) +
  labs(
    x = NULL, y = NULL, title = "Chromosome-level diversity statistics",
    subtitle = "Sim. points represent one chrom. mean value; empirical points represent one chrom. value"
    ) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
  theme_bw(base_size = 24) +
  theme(
    legend.position = "top", legend.direction = "horizontal",
    legend.title = element_blank(), panel.grid.minor = element_blank(),
    strip.background = element_blank(), strip.text = element_text(face = "bold")
    )
chrom.div.box.plot

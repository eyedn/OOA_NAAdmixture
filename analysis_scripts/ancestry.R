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
library(nanoparquet)
library(glue)
library(ggridges)
library(ggrepel)
library(broom)
library(lme4)
library(lmerTest)
library(emmeans)
library(boot)

set.seed(123)

sim.small.data.dir <- "~/scratch/OOA_NAAdmixture_small/stats"
sim.large.data.dir <- "~/scratch/OOA_NAAdmixture_large/stats"
emp.data.dir <- "~/scratch/OOA_NAAdmixture_1kG_cp/stats"

chrom.levels <- c(as.character(22:1), "all")
chrom.lengths <- read_table(
  "~/proj/1000GenomeNYGC_hg38_karatas/ONEKG_chr_lens.tsv", col_names = TRUE, 
  show_col_types = FALSE
  ) %>%
  mutate(chr = str_replace(as.character(chr), "chr", "")) %>%
  rename(chrom = chr)


# read in ancestry data ----


## small simulation ----


### small simulation: tspop ---- 
sim.small.anc.chr <- map_dfr(
  1:22,
  \(chr) {
    file.path(sim.small.data.dir, glue("ancestry.chr{chr}.parquet")) %>%
      read_parquet() %>%
      mutate(
        rep = as.numeric(rep), component_1_q = as.numeric(afr_tspop),
        component_2_q = as.numeric(eur_tspop), component_3_q = NA,
        component_4_q = NA, component_5_q = NA, span = as.numeric(span), k = 0,
        data.type = "Simulation_small", method = "tspop"
        ) %>%
      select(-afr_tspop, -eur_tspop, -afr_q, -eur_q)
    }
  )

# sim.small.anc.genome <- read_parquet(
#   file.path(sim.small.data.dir, glue("ancestry.parquet"))
#   ) %>%
#   mutate(
#     rep = as.numeric(rep), chrom = "all", component_1_q = as.numeric(afr_tspop),
#     component_2_q = as.numeric(eur_tspop), span = as.numeric(span),
#     component_3_q = NA, component_4_q = NA, component_5_q = NA, k = 0, 
#     data.type = "Simulation_small", method = "Truth"
#     ) %>%
#   select(-afr_tspop, -eur_tspop, -afr_q, -eur_q)


### small simulation: ADMIXTURE ---- 
sim.small.ADMIXTURE.anc.chr <- map_dfr(
  1:22,
  \(chr) {
    file.path(
      sim.small.data.dir, glue("ancestry_ADMIXTURE_multik.chr{chr}.parquet")
    ) %>%
      read_parquet() %>%
      mutate(
        chrom = as.character(chrom), sample_id = as.character(sample_id),
        data.type = "Simulation_small", method = "ADMIXTURE"
      )
    }
  )

# sim.small.ADMIXTURE.anc.genome <- read_parquet(
#   file.path(sim.small.data.dir, glue("ancestry_ADMIXTURE_multik.parquet"))
#   ) %>%
#   mutate(
#     chrom = "all", sample_id = as.character(sample_id), 
#     data.type = "Simulation_small", method = "ADMIXTURE"
#     )

### small simulation: fastStructure ---- 
sim.small.fastStructure.anc.chr <- map_dfr(
  1:22,
  \(chr) {
    file.path(
      sim.small.data.dir, glue("ancestry_fastStructure_multik.chr{chr}.parquet")
    ) %>%
      read_parquet() %>%
      mutate(
        chrom = as.character(chrom), sample_id = as.character(sample_id),
        data.type = "Simulation_small", method = "fastStructure"
      )
    }
  )

# sim.small.fastStructure.anc.genome <- read_parquet(
#   file.path(
#     sim.small.data.dir, glue("ancestry_fastStructure_multik.parquet")
#     )
#   ) %>%
#   mutate(
#     chrom = "all", sample_id = as.character(sample_id), 
#     data.type = "Simulation_small", method = "fastStructure"
#     )


## large simulation ----


### large simulation: tspop ---- 
sim.large.anc.chr <- map_dfr(
  1:22,
  \(chr) {
    file.path(sim.large.data.dir, glue("ancestry.chr{chr}.parquet")) %>%
      read_parquet() %>%
      mutate(
        rep = as.numeric(rep), component_1_q = as.numeric(afr_tspop),
        component_2_q = as.numeric(eur_tspop), component_3_q = NA,
        component_4_q = NA, component_5_q = NA, span = as.numeric(span), k = 0,
        data.type = "Simulation_large", method = "tspop"
      ) %>%
      select(-afr_tspop, -eur_tspop, -afr_q, -eur_q)
  }
)

# sim.large.anc.genome <- read_parquet(
#   file.path(sim.large.data.dir, glue("ancestry.parquet"))
#   ) %>%
#   mutate(
#     rep = as.numeric(rep), chrom = "all", component_1_q = as.numeric(afr_tspop),
#     component_2_q = as.numeric(eur_tspop), span = as.numeric(span),
#     component_3_q = NA, component_4_q = NA, component_5_q = NA, k = 0, 
#     data.type = "Simulation_large", method = "Truth"
#     ) %>%
#   select(-afr_tspop, -eur_tspop, -afr_q, -eur_q)


### large simulation: ADMIXTURE ---- 
sim.large.ADMIXTURE.anc.chr <- map_dfr(
  1:22,
  \(chr) {
    file.path(
      sim.large.data.dir, glue("ancestry_ADMIXTURE_multik.chr{chr}.parquet")
    ) %>%
      read_parquet() %>%
      mutate(
        chrom = as.character(chrom), sample_id = as.character(sample_id),
        data.type = "Simulation_large", method = "ADMIXTURE"
      )
  }
)

# sim.ADMIXTURE.anc.genome <- read_parquet(
#   file.path(sim.small.data.dir, glue("ancestry_ADMIXTURE_multik.parquet"))
#   ) %>%
#   mutate(
#     chrom = "all", sample_id = as.character(sample_id), 
#     data.type = "Simulation_small", method = "ADMIXTURE"
#     )

## small simulation: fastStructure ---- 
sim.large.fastStructure.anc.chr <- map_dfr(
  1:22,
  \(chr) {
    file.path(
      sim.large.data.dir, glue("ancestry_fastStructure_multik.chr{chr}.parquet")
    ) %>%
      read_parquet() %>%
      mutate(
        chrom = as.character(chrom), sample_id = as.character(sample_id),
        data.type = "Simulation_large", method = "fastStructure"
      )
  }
)

# sim.large.fastStructure.anc.genome <- read_parquet(
#   file.path(
#     sim.large.data.dir, glue("ancestry_fastStructure_multik.parquet")
#     )
#   ) %>%
#   mutate(
#     chrom = "all", sample_id = as.character(sample_id), 
#     data.type = "Simulation_large", method = "fastStructure"
#     )


## empirical ----


### empirical: ADMIXTURE ----
emp.ADMIXTURE.anc.chr <- map_dfr(
  1:22,
  \(chr) {
    file.path(
      emp.data.dir, glue("ancestry_ADMIXTURE_multik.chr{chr}.parquet")
      ) %>%
      read_parquet() %>%
      mutate(
        chrom = as.character(chrom), sample_id = as.character(sample_id),
        data.type = "Empirical", method = "ADMIXTURE"
      )
    }
  )

emp.ADMIXTURE.anc.genome <- read_parquet(
  file.path(emp.data.dir, glue("ancestry_ADMIXTURE_multik.parquet"))
  ) %>%
  mutate(
    chrom = "all", sample_id = as.character(sample_id), data.type = "Empirical", 
    method = "ADMIXTURE"
    )


### empirical: fastStructure ----
emp.fastStructure.anc.chr <- map_dfr(
  1:22,
  \(chr) {
    file.path(
      emp.data.dir, glue("ancestry_fastStructure_multik.chr{chr}.parquet")
      ) %>%
      read_parquet() %>%
      mutate(
        chrom = as.character(chrom), sample_id = as.character(sample_id),
        data.type = "Empirical", method = "fastStructure"
      )
    }
  )

emp.fastStructure.anc.genome <- read_parquet(
  file.path(
    emp.data.dir, glue("ancestry_fastStructure_multik.parquet")
    )
  ) %>%
  mutate(
    chrom = "all", sample_id = as.character(sample_id), data.type = "Empirical", 
    method = "fastStructure"
    )


## empirical: FLARE ----
# NOTE: results from traceAdmix project, saved as FLARE_anc_lengths.rds
# emp.FLARE.anc.chr <- readRDS("~/proj/traceAdmix/FLARE_anc_lengths.rds") %>%
#   rename(chrom = chr, vcf_sample_id = sam, component_1_q = anc1.frac) %>%
#   filter(vcf_sample_id %in% emp.ADMIXTURE.anc.chr$vcf_sample_id) %>%
#   mutate(
#     chrom = as.character(chrom), rep = 0, sample_id = NA, k = 0,
#     component_2_q = 1 - component_1_q, component_3_q = NA, component_4_q = NA,
#     component_5_q = NA, span = NA, data.type = "Empirical", method = "FLARE"
#     ) %>%
#   select(colnames(emp.ADMIXTURE.anc.chr))


## combine ancestry data ----
anc <- rbind(
  sim.small.anc.chr, sim.small.ADMIXTURE.anc.chr, 
  sim.small.fastStructure.anc.chr, sim.large.anc.chr, 
  sim.large.ADMIXTURE.anc.chr, sim.large.fastStructure.anc.chr,
  emp.ADMIXTURE.anc.chr, emp.fastStructure.anc.chr, emp.ADMIXTURE.anc.genome, 
  emp.fastStructure.anc.genome
  ) %>%
  mutate(
    role = ifelse(
      pop %in% c("AFR", "YRI"), "AFR",
      ifelse(pop %in% c("EUR", "CEU"), "EUR", "ADX")
      )
    )


# determine ancestry components ----


## inspect K == 2 for simulated ancestry ----
sim.admixPlot.df <- anc %>%
  filter(pop == "ADX") %>%
  mutate(chrom = factor(chrom, levels = chrom.levels)) %>%
  pivot_longer(
    starts_with("component_"), names_to = "component", values_to = "q"
    ) %>%
  filter(component %in% c("component_1_q", "component_2_q")) 

sim.admixPlot.summary <- sim.admixPlot.df %>%
  group_by(chrom, rep, data.type, method, component) %>%
  summarise(
    mean_q = mean(q, na.rm = TRUE), sd_q = sd(q, na.rm = TRUE), .groups = "drop"
    )
sim.admixPlot.summary


### small simulation: tspop ----
sim.small.rep1.tspop.admixPlot <- ggplot(
  sim.admixPlot.df %>% 
    filter(method == "tspop", data.type == "Simulation_small", rep == 1),
  aes(x = vcf_sample_id, y = q, fill = component)
  ) +
  geom_col() +
  facet_wrap(~chrom, ncol = 4) +
  labs(
    x = NULL, y = "Ancestry proportion",
    title = "fastStructure Admixture plots"
  ) +
  theme_bw(base_size = 24) +
  theme(
    axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    strip.background = element_blank(), panel.spacing.x = unit(0.1, "lines")
  )
sim.small.rep1.tspop.admixPlot


### small simulation: ADMIXTURE ----
sim.small.rep1.ADMIXTURE.admixPlot <- ggplot(
  sim.admixPlot.df %>% 
    filter(method == "ADMIXTURE", rep == 1, data.type == "Simulation_small"),
  aes(x = vcf_sample_id, y = q, fill = component)
  ) +
  geom_col() +
  facet_wrap(~chrom, ncol = 4, scales = "free_x") +
  labs(
    x = NULL, y = "Ancestry proportion",
    title = "ADMIXTURE Admixture plots"
  ) +
  theme_bw(base_size = 24) +
  theme(
    axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    strip.background = element_blank(), panel.spacing.x = unit(0.1, "lines")
  )
sim.small.rep1.ADMIXTURE.admixPlot


### small simulation: fastStructure ----
sim.small.rep1.fastStructure.admixPlot <- ggplot(
  sim.admixPlot.df %>% 
    filter(method == "fastStructure", rep == 1, data.type == "Simulation_small"),
  aes(x = vcf_sample_id, y = q, fill = component)
) +
  geom_col() +
  facet_wrap(~chrom, ncol = 4, scales = "free_x") +
  labs(
    x = NULL, y = "Ancestry proportion",
    title = "fastStructure Admixture plots"
  ) +
  theme_bw(base_size = 24) +
  theme(
    axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    strip.background = element_blank(), panel.spacing.x = unit(0.1, "lines")
  )
sim.small.rep1.fastStructure.admixPlot


### small simulation: tspop ----
sim.large.rep1.tspop.admixPlot <- ggplot(
  sim.admixPlot.df %>% 
    filter(method == "tspop", rep == 1, data.type == "Simulation_large"),
  aes(x = vcf_sample_id, y = q, fill = component)
) +
  geom_col() +
  facet_wrap(~chrom, ncol = 4) +
  labs(
    x = NULL, y = "Ancestry proportion",
    title = "fastStructure Admixture plots"
  ) +
  theme_bw(base_size = 24) +
  theme(
    axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    strip.background = element_blank(), panel.spacing.x = unit(0.1, "lines")
  )
sim.large.rep1.tspop.admixPlot


### large simulation: ADMIXTURE ----
sim.large.rep1.ADMIXTURE.admixPlot <- ggplot(
  sim.admixPlot.df %>% 
    filter(method == "ADMIXTURE", rep == 1, data.type == "Simulation_large"),
  aes(x = vcf_sample_id, y = q, fill = component)
) +
  geom_col() +
  facet_wrap(~chrom, ncol = 4, scales = "free_x") +
  labs(
    x = NULL, y = "Ancestry proportion",
    title = "ADMIXTURE Admixture plots"
  ) +
  theme_bw(base_size = 24) +
  theme(
    axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    strip.background = element_blank(), panel.spacing.x = unit(0.1, "lines")
  )
sim.large.rep1.ADMIXTURE.admixPlot


### large simulation: fastStructure ----
sim.large.rep1.fastStructure.admixPlot <- ggplot(
  sim.admixPlot.df %>% 
    filter(method == "fastStructure", rep == 1, data.type == "Simulation_large"),
  aes(x = vcf_sample_id, y = q, fill = component)
) +
  geom_col() +
  facet_wrap(~chrom, ncol = 4, scales = "free_x") +
  labs(
    x = NULL, y = "Ancestry proportion",
    title = "fastStructure Admixture plots"
  ) +
  theme_bw(base_size = 24) +
  theme(
    axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    strip.background = element_blank(), panel.spacing.x = unit(0.1, "lines")
  )
sim.large.rep1.fastStructure.admixPlot


## inspect optimal k for empirical fastStructure and ADMIXTURE ----


### read in fastStructure chooseK results for empirical data ----
emp.fastStructure.choseK.chr <- map_dfr(
  1:22,
  \(chr) {
    file.path(
      emp.data.dir, glue("fastStructure_chooseK.chr{chr}.parquet")
    ) %>%
      read_parquet()
  }
)
table(emp.fastStructure.choseK.chr$max_marginal_likelihood_k)
table(emp.fastStructure.choseK.chr$model_components_k)

emp.fastStructure.choseK.genom <- read_parquet(
  file.path(emp.data.dir, glue("fastStructure_chooseK.parquet"))
)
table(emp.fastStructure.choseK.genom$max_marginal_likelihood_k)
table(emp.fastStructure.choseK.genom$model_components_k)

### assign optimal k ----
k.use <- 2
emp.admixPlot.df <- anc %>%
  filter(pop == "ASW", k %in% c(0, k.use)) %>%
  mutate(chrom = factor(chrom, levels = chrom.levels)) %>%
  pivot_longer(
    starts_with("component_"), names_to = "component", values_to = "q"
  ) %>%
  filter(component %in% paste0("component_", seq_len(k.use), "_q"))

emp.admixPlot.summary <- emp.admixPlot.df %>%
  group_by(chrom, method, component) %>%
  summarise(
    mean_q = mean(q, na.rm = TRUE), sd_q = sd(q, na.rm = TRUE), .groups = "drop"
    )
emp.admixPlot.summary


### empirical: ADMIXTURE ----
emp.ADMIXTURE.admixPlot <- ggplot(
  emp.admixPlot.df %>% filter(method == "ADMIXTURE"),
  aes(x = vcf_sample_id, y = q, fill = component)
) +
  geom_col() +
  facet_wrap(~chrom, ncol = 4, scales = "free_x") +
  labs(
    x = NULL, y = "Ancestry proportion",
    title = "ADMIXTURE Admixture plots"
  ) +
  theme_bw(base_size = 24) +
  theme(
    axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    strip.background = element_blank(), panel.spacing.x = unit(0.1, "lines")
  )
emp.ADMIXTURE.admixPlot


### empirical: fastStructure
emp.fastStructure.admixPlot <- ggplot(
  emp.admixPlot.df %>% filter(method == "fastStructure"),
  aes(x = vcf_sample_id, y = q, fill = component)
) +
  geom_col() +
  facet_wrap(~chrom, ncol = 4, scales = "free_x") +
  labs(
    x = NULL, y = "Ancestry proportion",
    title = "fastStructure Admixture plots"
  ) +
  theme_bw(base_size = 24) +
  theme(
    axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    strip.background = element_blank(), panel.spacing.x = unit(0.1, "lines")
  )
emp.fastStructure.admixPlot


### empirical: FLARE
# emp.FLARE.admixPlot <- ggplot(
#   emp.admixPlot.df %>% filter(method == "FLARE"),
#   aes(x = vcf_sample_id, y = q, fill = component)
# ) +
#   geom_col() +
#   facet_wrap(~chrom, ncol = 4, scales = "free_x") +
#   labs(
#     x = NULL, y = "Ancestry proportion",
#     title = "fastStructure Admixture plots"
#   ) +
#   theme_bw(base_size = 24) +
#   theme(
#     axis.text.x = element_blank(), axis.ticks.x = element_blank(),
#     strip.background = element_blank(), panel.spacing.x = unit(0.1, "lines")
#   )
# emp.FLARE.admixPlot


## determine which component is AFR based on max average ancestry ----
component.map <- anc %>%
  filter(role == "ADX", k %in% c(0, k.use)) %>%
  group_by(rep, pop, chrom, data.type, method) %>%
  summarise(
    mean.1 = mean(component_1_q, na.rm = TRUE),
    sd.1 = sd(component_1_q, na.rm = TRUE),
    mean.2 = mean(component_2_q, na.rm = TRUE),
    sd.2 = sd(component_2_q, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    afr.component = if_else(
      mean.1 >= mean.2,
      "component_1_q", "component_2_q"
    )
  )

anc <- anc %>%
  filter(k %in% c(0, k.use)) %>%
  left_join(
    component.map %>% select(rep, chrom, data.type, method, afr.component),
    by = c("rep", "chrom", "data.type", "method")
  ) %>%
  mutate(
    afr.q = case_when(
      afr.component == "component_1_q" ~ component_1_q,
      afr.component == "component_2_q" ~ component_2_q,
      TRUE ~ NA_real_
    ),
    eur.q = case_when(
      afr.component == "component_1_q" ~ component_2_q,
      afr.component == "component_2_q" ~ component_1_q,
      TRUE ~ NA_real_
    )
  ) %>%
  select(
    rep, chrom, pop, sample_id, vcf_sample_id, span, data.type, method, role, 
    afr.q, eur.q
    )


# summarize data ----
boot.mean <- function(data, indices) {
  x <- data[indices]
  return(mean = mean(x))
  }

boot.sd <- function(data, indices) {
  x <- data[indices]
  return(sd = sd(x))
  }

anc.summary <- anc %>%
  filter(role == "ADX") %>%
  group_by(rep, chrom, pop, data.type, method) %>%
  group_modify(
    ~ {
      x <- .x$afr.q
      x <- x[!is.na(x)]
      
      # default for groups that will not be bootstrapped
      mean.boot       <- NA_real_
      mean.boot.lower <- NA_real_
      mean.boot.upper <- NA_real_
      sd.boot         <- NA_real_
      sd.boot.lower   <- NA_real_
      sd.boot.upper   <- NA_real_
      
      # bootstrap only empirical groups
      if (.y$data.type[[1]] == "Empirical") {
        boot.mean.fit <- boot(data = x, statistic = boot.mean, R = 1000)
        boot.mean.ci <- boot.ci(boot.mean.fit, type = "bca")
        boot.sd.fit <- boot(data = x, statistic = boot.sd, R = 1000)
        boot.sd.ci <- boot.ci(boot.sd.fit, type = "bca")
        
        mean.boot       <- boot.mean.fit$t0
        mean.boot.lower <- boot.mean.ci$bca[4]
        mean.boot.upper <- boot.mean.ci$bca[5]
        sd.boot         <- boot.sd.fit$t0
        sd.boot.lower   <- boot.sd.ci$bca[4]
        sd.boot.upper   <- boot.sd.ci$bca[5]
      }
      
      tibble(
        mean = mean(x),
        mean.boot = mean.boot,
        mean.boot.lower = mean.boot.lower,
        mean.boot.upper = mean.boot.upper,
        sd = sd(x),
        sd.boot = sd.boot,
        sd.boot.lower = sd.boot.lower,
        sd.boot.upper = sd.boot.upper,
        median = median(x),
        q25 = quantile(x, probs = 0.25, names = FALSE),
        q75 = quantile(x, probs = 0.75, names = FALSE),
        n = length(x)
        )
      }
    ) %>%
  ungroup()
anc.summary

n.downsample <- 50
n.iterations <- 1000
anc.sim.downsample.summary <- anc %>%
  filter(role == "ADX", data.type != "Empirical") %>%
  group_by(rep, chrom, pop, data.type, method) %>%
  group_modify(
    ~ {
      x <- .x$afr.q
      x <- x[!is.na(x)]
      x.sample <- sample(x, size = n.downsample, replace = FALSE)

      tibble(
        # average estimate across downsampling iterations
        mean.rand.downsample = mean(x.sample),
        sd.rand.downsample = sd(x.sample),
        median.rand.downsample = median(x.sample),
        q25.rand.downsample = quantile(x.sample, probs = 0.25, names = FALSE),
        q75.rand.downsample = quantile(x.sample, probs = 0.75, names = FALSE),
        
        n.rand.downsample = n.downsample,
        n.iterations = n.iterations
      )
    }
  ) %>%
  ungroup()
anc.sim.downsample.summary

anc.summary <- anc.summary %>%
  left_join(
    anc.sim.downsample.summary,
    by = c("rep", "chrom", "pop", "data.type", "method")
    )
anc.summary

anc.summary.chrom.levels <- c("1", "5", "10", "14", "18", "22")
anc.summary.plot.df <- anc.summary %>%
  filter(
    (
      chrom %in% anc.summary.chrom.levels &
        data.type == "Simulation_small"
      ) |
      (
        chrom %in% anc.summary.chrom.levels &
          data.type == "Simulation_large"
        ) |
      (
        chrom %in% c("all", anc.summary.chrom.levels) &
         data.type == "Empirical" & method != "FLARE"
        )
    ) %>%
  mutate(
    chrom = factor(
      chrom, levels = c("all", anc.summary.chrom.levels)
      ),
    method = factor(
      method,
      levels = c("tspop", "ADMIXTURE", "fastStructure")
    )
  )

## boxplots of mean ancestry ----
make.mean.ancestry.plot <- function(
    plot.method, empirical.color, chromosome.levels
    ) {
  
  # simulation boxplots: truth plus selected inference method
  simulation.df <- anc.summary.plot.df %>%
    filter(
      data.type != "Empirical", method %in% c("tspop", plot.method),
      chrom %in% chromosome.levels
      )
  
  # chromosome-specific empirical estimates
  empirical.chrom.df <- anc.summary.plot.df %>%
    filter(
      data.type == "Empirical", method == plot.method,
      chrom %in% chromosome.levels
      )
  
  # whole-genome empirical estimate
  empirical.genome.df <- anc.summary.plot.df %>%
    filter(
      data.type == "Empirical", method == plot.method, chrom == "all"
      )
  
  fill.values <- c(
    "tspop.Simulation_small" = "#00AEDB", 
    "tspop.Simulation_large" = "#007E9F",
    "ADMIXTURE.Simulation_small" = "#FA9EBC",
    "ADMIXTURE.Simulation_large" = "#E66791",
    "fastStructure.Simulation_small" = "#FFDBD1",
    "fastStructure.Simulation_large" = "#E69B88"
    )
  
  fill.breaks <- if (plot.method == "ADMIXTURE") {
    c(
      "tspop.Simulation_small", "tspop.Simulation_large",
      "ADMIXTURE.Simulation_small", "ADMIXTURE.Simulation_large"
      )
    } else {
    c(
      "tspop.Simulation_small", "tspop.Simulation_large",
      "fastStructure.Simulation_small", "fastStructure.Simulation_large"
      )
    }
  
  fill.labels <- if (plot.method == "ADMIXTURE") {
    c(
      "Truth (small)", "Truth (large)", "ADMIXTURE (small)", "ADMIXTURE (large)"
      )
    } else {
    c(
      "Truth (small)", "Truth (large)", "fastStructure (small)",
      "fastStructure (large)"
      )
    }
  
  ggplot(
    simulation.df,
    aes(x = chrom, y = mean, fill = interaction(method, data.type))
    ) +
    # whole-genome empirical 95% CI
    geom_rect(
      data = empirical.genome.df,
      aes(
        xmin = -Inf, xmax = Inf, ymin = mean.boot.lower, ymax = mean.boot.upper
        ),
      inherit.aes = FALSE, fill = empirical.color, color = empirical.color,
      alpha = 0.10, linetype = "dotted"
    ) +
    # whole-genome empirical mean
    geom_hline(
      data = empirical.genome.df, aes(yintercept = mean.boot),
      inherit.aes = FALSE, color = empirical.color, linetype = "dashed", 
      linewidth = 0.8
      ) +
    # simulation values
    geom_boxplot(
      aes(group = interaction(chrom, data.type, method)), outlier.shape = NA,
      color = "black"
      ) +
    # chromosome-specific empirical 95% CI
    # geom_errorbar(
    #   data = empirical.chrom.df,
    #   aes(x = chrom, ymin = mean.boot.lower, ymax = mean.boot.upper),
    #   inherit.aes = FALSE, color = empirical.color, width = 0.15, 
    #   linewidth = 1.1
    #   ) +
    # chromosome-specific empirical mean
    # geom_point(
    #   data = empirical.chrom.df,
    #   aes(x = chrom, y = mean.boot), inherit.aes = FALSE, shape = 23, size = 4,
    #   stroke = 1.1, color = "black", fill = empirical.color
    # ) +
    # formatting
    scale_fill_manual(
      name = NULL, breaks = fill.breaks, labels = fill.labels, 
      values = fill.values
      ) +
    scale_x_discrete(
      limits = chromosome.levels,
      drop = FALSE
    ) +
    labs(
      x = "Chromosome",
      y = "Mean African ancestry",
      title = paste(
        "Sample Mean of African Ancestry:",
        plot.method
      ),
      subtitle = paste(
        "Boxplots: simulation means;",
        "diamonds and error bars: empirical chromosome means and 95% CIs;",
        "dashed line and band: empirical whole-genome mean and 95% CI", 
        sep = "\n"
      )
    ) +
    guides(
      fill = guide_legend(
        order = 1,
        override.aes = list(
          shape = 23, color = "black", linewidth = 1, size = 3
          )
        )
      ) +
    theme_bw(base_size = 24) +
    theme(legend.position = "top", panel.grid.minor = element_blank())
  }

meanAfrAnc.admixture.plot <- make.mean.ancestry.plot(
  plot.method = "ADMIXTURE", empirical.color = "#B83264", 
  chromosome.levels = anc.summary.chrom.levels
  )
meanAfrAnc.admixture.plot

meanAfrAnc.fastStructure.plot <- make.mean.ancestry.plot(
  plot.method = "fastStructure", empirical.color = "#B9584A", 
  chromosome.levels = anc.summary.chrom.levels
  )
meanAfrAnc.fastStructure.plot


## boxplots of sd ancestry ----
make.sd.ancestry.plot <- function(
    plot.method, empirical.color, chromosome.levels
    ) {
  
  # truth plus the selected inference method
  simulation.df <- anc.summary.plot.df %>%
    filter(
      data.type != "Empirical", method %in% c("tspop", plot.method),
      chrom %in% chromosome.levels
      )
  
  # chromosome-specific empirical estimates
  empirical.chrom.df <- anc.summary.plot.df %>%
    filter(
      data.type == "Empirical", method == plot.method,
      chrom %in% chromosome.levels
      )
  
  # Whole-genome empirical estimate
  empirical.genome.df <- anc.summary.plot.df %>%
    filter(
      data.type == "Empirical", method == plot.method, chrom == "all"
      )
  
  fill.values <- c(
    "tspop.Simulation_small" = "#00AEDB",
    "tspop.Simulation_large" = "#007E9F",
    "ADMIXTURE.Simulation_small" = "#FA9EBC",
    "ADMIXTURE.Simulation_large" = "#E66791",
    "fastStructure.Simulation_small" = "#FFDBD1",
    "fastStructure.Simulation_large" = "#E69B88"
  )
  
  if (plot.method == "ADMIXTURE") {
    fill.breaks <- c(
      "tspop.Simulation_small", "tspop.Simulation_large",
      "ADMIXTURE.Simulation_small", "ADMIXTURE.Simulation_large"
    )
    
    fill.labels <- c(
      "Truth (small)", "Truth (large)", "ADMIXTURE (small)", "ADMIXTURE (large)"
      )
    } else {
    fill.breaks <- c(
      "tspop.Simulation_small", "tspop.Simulation_large", 
      "fastStructure.Simulation_small", "fastStructure.Simulation_large"
      )
    
    fill.labels <- c(
      "Truth (small)", "Truth (large)", "fastStructure (small)", 
      "fastStructure (large)"
      )
    }
  
  ggplot(
    simulation.df,
    aes(x = chrom, y = sd, fill = interaction(method, data.type))
    ) +
    # whole-genome empirical 95% CI
    geom_rect(
      data = empirical.genome.df,
      aes(xmin = -Inf, xmax = Inf, ymin = sd.boot.lower, ymax = sd.boot.upper),
      inherit.aes = FALSE, fill = empirical.color, color = empirical.color,
      alpha = 0.10, linetype = "dotted"
      ) +
    # whole-genome empirical SD
    geom_hline(
      data = empirical.genome.df, aes(yintercept = sd.boot), 
      inherit.aes = FALSE, color = empirical.color, linetype = "dashed",
      linewidth = 0.8
      ) +
    # simulation SD values
    geom_boxplot(
      aes(group = interaction(chrom, data.type, method)), outlier.shape = NA, 
      color = "black"
      ) +
    # chromosome-specific empirical 95% CIs
    geom_errorbar(
      data = empirical.chrom.df,
      aes(x = chrom, ymin = sd.boot.lower, ymax = sd.boot.upper),
      inherit.aes = FALSE, color = empirical.color, width = 0.15, 
      linewidth = 1.1
    ) +
    # chromosome-specific empirical SDs
    geom_point(
      data = empirical.chrom.df, aes(x = chrom, y = sd.boot), 
      inherit.aes = FALSE, shape = 23, size = 4, stroke = 1.1,
      color = "black", fill = empirical.color
    ) +
    # formatting
    scale_fill_manual(
      name = NULL, breaks = fill.breaks, labels = fill.labels, 
      values = fill.values
      ) +
    scale_x_discrete(limits = chromosome.levels, drop = FALSE) +
    labs(
      x = "Chromosome",
      y = "Standard deviation of African ancestry",
      title = paste(
        "Sample Standard Deviation of African Ancestry:",
        plot.method
      ),
      subtitle = paste(
        "Boxplots: simulation standard deviations;",
        "diamonds and error bars: empirical chromosome SDs and 95% CIs;",
        "dashed line and band: empirical whole-genome SD and 95% CI",
        sep = "\n"
      )
    ) +
    guides(
      fill = guide_legend(
        order = 1,
        override.aes = list(
          shape = 23, color = "black", linewidth = 1, size = 3
          )
        )
      ) +
    theme_bw(base_size = 24) +
    theme(legend.position = "top", panel.grid.minor = element_blank())
  }


sdAfrAnc.admixture.plot <- make.sd.ancestry.plot(
  plot.method = "ADMIXTURE", empirical.color = "#B83264", 
  chromosome.levels = anc.summary.chrom.levels
  )
sdAfrAnc.admixture.plot

sdAfrAnc.fastStructure.plot <- make.sd.ancestry.plot(
  plot.method = "fastStructure", empirical.color = "#B9584A",
  chromosome.levels = anc.summary.chrom.levels
  )
sdAfrAnc.fastStructure.plot


## combined boxplot for mean and sd ----
anc.summary.facet.df <- anc.summary.plot.df %>%
  filter(data.type != "Empirical") %>%
  pivot_longer(
    cols = c(mean, sd), names_to = "stat", values_to = "estimate"
    ) %>%
  mutate(
    stat = factor(
      stat, levels = c("mean", "sd"),
      labels = c("Mean ancestry", "Standard dev. of ancestry")
      )
    )
anc.empirical.lines.df <- anc.summary.plot.df %>%
  filter(
    data.type == "Empirical", chrom == "all", 
    method %in% c("ADMIXTURE", "fastStructure")
    ) %>%
  pivot_longer(
    cols = c(mean, sd), names_to = "stat", values_to = "estimate"
    ) %>%
  mutate(
    stat = factor(
      stat, levels = c("mean", "sd"),
      labels = c("Mean ancestry", "Standard dev. of ancestry"))
    )

# combined plot
make.afr.ancestry.facet.plot <- function(
    plot.method, empirical.color, chromosome.levels
    ) {
  
  stat.levels <- c("mean", "sd")
  stat.labels <- c(mean = "Mean ancestry", sd = "Standard dev. of ancestry")
  
  # simulation values: truth plus selected inference method
  simulation.facet.df <- anc.summary.plot.df %>%
    filter(
      data.type != "Empirical", method %in% c("tspop", plot.method),
      chrom %in% chromosome.levels
      ) %>%
    pivot_longer(
      cols = c(mean, sd), names_to = "stat", values_to = "estimate"
      ) %>%
    mutate(stat = factor(stat, levels = stat.levels))
  
  # chromosome-specific empirical estimates and 95% CIs
  empirical.chrom.facet.df <- bind_rows(
    anc.summary.plot.df %>%
      filter(
        data.type == "Empirical", method == plot.method,
        chrom %in% chromosome.levels
        ) %>%
      transmute(
        chrom, stat = "mean", estimate = mean.boot, lower = mean.boot.lower,
        upper = mean.boot.upper
        ),
    anc.summary.plot.df %>%
      filter(
        data.type == "Empirical", method == plot.method,
        chrom %in% chromosome.levels
        ) %>%
      transmute(
        chrom, stat = "sd", estimate = sd.boot, lower = sd.boot.lower,
        upper = sd.boot.upper
        )
    ) %>%
    mutate(stat = factor(stat, levels = stat.levels))
  
  # whole-genome empirical estimates and 95% CIs
  empirical.genome.facet.df <- bind_rows(
    anc.summary.plot.df %>%
      filter(
        data.type == "Empirical", method == plot.method, chrom == "all"
        ) %>%
      transmute(
        stat = "mean", estimate = mean.boot, lower = mean.boot.lower,
        upper = mean.boot.upper
        ),
    anc.summary.plot.df %>%
      filter(
        data.type == "Empirical", method == plot.method, chrom == "all"
        ) %>%
      transmute(
        stat = "sd", estimate = sd.boot, lower = sd.boot.lower,
        upper = sd.boot.upper
        )
    ) %>%
    mutate(stat = factor(stat, levels = stat.levels))
  
  fill.values <- c(
    "tspop.Simulation_small" = "#00AEDB",
    "tspop.Simulation_large" = "#007E9F",
    "ADMIXTURE.Simulation_small" = "#FA9EBC",
    "ADMIXTURE.Simulation_large" = "#E66791",
    "fastStructure.Simulation_small" = "#FFDBD1",
    "fastStructure.Simulation_large" = "#E69B88"
  )
  
  if (plot.method == "ADMIXTURE") {
    fill.breaks <- c(
      "tspop.Simulation_small", "tspop.Simulation_large",
      "ADMIXTURE.Simulation_small", "ADMIXTURE.Simulation_large"
      )
    fill.labels <- c(
      "Truth (small)", "Truth (large)", "ADMIXTURE (small)", "ADMIXTURE (large)"
      )
    } else {
    fill.breaks <- c(
      "tspop.Simulation_small", "tspop.Simulation_large",
      "fastStructure.Simulation_small", "fastStructure.Simulation_large"
      )
    fill.labels <- c(
      "Truth (small)", "Truth (large)", "fastStructure (small)", 
      "fastStructure (large)"
      )
    }
  
  ggplot(
    simulation.facet.df,
    aes(x = chrom, y = estimate, fill = interaction(method, data.type))
    ) +
    # whole-genome empirical 95% CI bands
    geom_rect(
      data = empirical.genome.facet.df,
      aes(xmin = -Inf, xmax = Inf, ymin = lower, ymax = upper),
      inherit.aes = FALSE, fill = empirical.color, color = empirical.color,
      alpha = 0.10, linetype = "dotted"
      ) +
    # whole-genome empirical estimates
    geom_hline(
      data = empirical.genome.facet.df, aes(yintercept = estimate),
      inherit.aes = FALSE, color = empirical.color, linetype = "dashed",
      linewidth = 0.8
      ) +
    # simulation estimates
    geom_boxplot(
      aes(group = interaction(chrom, data.type, method)), outlier.shape = NA,
      color = "black"
      ) +
    # chromosome-specific empirical 96% CIs
    geom_errorbar(
      data = empirical.chrom.facet.df,
      aes(x = chrom, ymin = lower, ymax = upper), inherit.aes = FALSE,
      color = empirical.color, width = 0.15, linewidth = 1.1
      ) +
    # chromosome-specific empirical estimates
    geom_point(
      data = empirical.chrom.facet.df, aes(x = chrom, y = estimate),
      inherit.aes = FALSE, shape = 23, size = 4, stroke = 1.1, color = "black",
      fill = empirical.color
      ) +
    # formating
    facet_grid(
      rows = vars(stat), scales = "free_y",
      labeller = labeller(stat = stat.labels)
      ) +
    scale_x_discrete(limits = chromosome.levels, drop = FALSE) +
    scale_fill_manual(
      name = NULL, breaks = fill.breaks, labels = fill.labels, 
      values = fill.values
      ) +
    labs(
      x = "Chromosome", y = NULL,
      title = paste(
        "African Ancestry on Admixed Chromosomes:",
        plot.method
      ),
      subtitle = paste(
        "Boxplots: simulation replicate estimates;",
        "diamonds and error bars: empirical chromosome estimates and 95% CIs;",
        "dashed lines and bands: empirical whole-genome estimates and 95% CIs",
        sep = "\n"
        )
      ) +
    guides(
      fill = guide_legend(
        order = 1,
        override.aes = list(shape = 23, color = "black", linewidth = 1,
          size = 3
          )
        )
      ) +
    theme_bw(base_size = 24) +
    theme(
      legend.position = "top", panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey95"),
      strip.text.y = element_text(size = 18)
    )
  }


afrAnc.admixture.facet.plot <- make.afr.ancestry.facet.plot(
  plot.method = "ADMIXTURE", empirical.color = "#B83264", 
  chromosome.levels = anc.summary.chrom.levels
  )
afrAnc.admixture.facet.plot

afrAnc.fastStructure.facet.plot <- make.afr.ancestry.facet.plot(
  plot.method = "fastStructure", empirical.color = "#B9584A",
  chromosome.levels = anc.summary.chrom.levels
  )
afrAnc.fastStructure.facet.plot


# mean ancestry vs chromosome length ----
anc.length.mean.df <- anc.summary %>%
  filter(chrom != "all") %>%
  left_join(chrom.lengths, by = "chrom") %>%
  mutate(
    chr.len.mb = if_else(
      data.type != "Empirical", chr_len, chr_len_after_qc
      ) / 1e6,
    ) %>%
  group_by(data.type) %>%
  mutate(chr.len.mb.c = chr.len.mb - mean(chr.len.mb)) %>%
  ungroup()


## plot ----
make.anc.length.mean.plot <- function(plot.method) {
  
  # Method-specific colors
  if (plot.method == "ADMIXTURE") {
    
    series.breaks <- c(
      "tspop.Simulation_small",
      "tspop.Simulation_large",
      "ADMIXTURE.Simulation_small",
      "ADMIXTURE.Simulation_large",
      "ADMIXTURE.Empirical"
    )
    
    series.labels <- c(
      "Truth (small)",
      "Truth (large)",
      "ADMIXTURE (small)",
      "ADMIXTURE (large)",
      "ADMIXTURE (empirical)"
    )
    
    series.colors <- c(
      "tspop.Simulation_small" = "#00AEDB",
      "tspop.Simulation_large" = "#007E9F",
      "ADMIXTURE.Simulation_small" = "#FA9EBC",
      "ADMIXTURE.Simulation_large" = "#E66791",
      "ADMIXTURE.Empirical" = "#B83264"
    )
    
    method.shape <- 25
    
  } else if (plot.method == "fastStructure") {
    
    series.breaks <- c(
      "tspop.Simulation_small",
      "tspop.Simulation_large",
      "fastStructure.Simulation_small",
      "fastStructure.Simulation_large",
      "fastStructure.Empirical"
    )
    
    series.labels <- c(
      "Truth (small)",
      "Truth (large)",
      "fastStructure (small)",
      "fastStructure (large)",
      "fastStructure (empirical)"
    )
    
    series.colors <- c(
      "tspop.Simulation_small" = "#00AEDB",
      "tspop.Simulation_large" = "#007E9F",
      "fastStructure.Simulation_small" = "#FFDBD1",
      "fastStructure.Simulation_large" = "#E69B88",
      "fastStructure.Empirical" = "#B9584A"
    )
    
    method.shape <- 24
    
  } else {
    stop("plot.method must be 'ADMIXTURE' or 'fastStructure'")
  }
  
  # Truth plus the selected inference method
  plot.df <- anc.length.mean.df %>%
    filter(
      (
        data.type %in% c("Simulation_small", "Simulation_large") &
          method %in% c("tspop", plot.method)
      ) |
        (
          data.type == "Empirical" &
            method == plot.method
        )
    ) %>%
    group_by(
      chrom,
      pop,
      data.type,
      method,
      chr.len.mb
    ) %>%
    summarize(
      estimate = median(mean, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      chrom = as.character(chrom),
      method = factor(
        method,
        levels = c("tspop", "ADMIXTURE", "fastStructure")
      ),
      data.type = factor(
        data.type,
        levels = c(
          "Simulation_small",
          "Simulation_large",
          "Empirical"
        )
      ),
      series = interaction(
        method,
        data.type,
        sep = ".",
        drop = TRUE
      ),
      point.shape = if_else(
        method == "tspop",
        21,
        method.shape
      ),
      label = (
        data.type == "Empirical" &
          as.integer(chrom) %% 3 == 1
      )
    )
  
  ggplot(
    plot.df,
    aes(
      x = chr.len.mb,
      y = estimate,
      fill = series,
      color = series,
      linetype = data.type,
      group = series
    )
  ) +
    geom_smooth(
      method = "lm",
      formula = y ~ x,
      se = FALSE,
      linewidth = 2
    ) +
    geom_point(
      aes(shape = method),
      color = "black",
      stroke = 1,
      size = 4
    ) +
    geom_text_repel(
      data = plot.df %>% filter(label),
      aes(label = chrom),
      inherit.aes = FALSE,
      x = plot.df %>% filter(label) %>% pull(chr.len.mb),
      y = plot.df %>% filter(label) %>% pull(estimate),
      color = "black",
      size = 5,
      seed = 123,
      min.segment.length = 0,
      force = 15,
      force_pull = 5,
      box.padding = 0.5,
      point.padding = 0.25,
      max.iter = 10000,
      show.legend = FALSE
    ) +
    scale_fill_manual(
      name = NULL,
      breaks = series.breaks,
      labels = series.labels,
      values = series.colors
    ) +
    scale_color_manual(
      name = NULL,
      breaks = series.breaks,
      labels = series.labels,
      values = series.colors
    ) +
    scale_shape_manual(
      name = NULL,
      values = c(
        "tspop" = 21,
        "ADMIXTURE" = 25,
        "fastStructure" = 24
      )
    ) +
    scale_linetype_manual(
      name = NULL,
      values = c(
        "Simulation_small" = "solid",
        "Simulation_large" = "longdash",
        "Empirical" = "dashed"
      )
    ) +
    labs(
      x = "Chromosome length (Mb)",
      y = "Mean African ancestry",
      title = paste(
        "Chromosome Ancestry versus Chromosome Length:",
        plot.method
      )
    ) +
    guides(
      color = "none",
      shape = "none",
      linetype = "none",
      fill = guide_legend(
        order = 1,
        override.aes = list(
          color = "black",
          linetype = 0,
          shape = c(
            21,
            21,
            method.shape,
            method.shape,
            method.shape
          )
        )
      )
    ) +
    theme_bw(base_size = 24) +
    theme(
      legend.position = "top",
      panel.grid.minor = element_blank()
    )
}

anc.length.mean.admixture.plot <- make.anc.length.mean.plot(
  plot.method = "ADMIXTURE"
)
anc.length.mean.admixture.plot

anc.length.mean.fastStructure.plot <- make.anc.length.mean.plot(
  plot.method = "fastStructure"
)
anc.length.mean.fastStructure.plot


# sd ancestry vs chromosome length ----
make.anc.length.sd.plot <- function(plot.method) {
  
  if (plot.method == "ADMIXTURE") {
    
    series.breaks <- c(
      "tspop.Simulation_small",
      "tspop.Simulation_large",
      "ADMIXTURE.Simulation_small",
      "ADMIXTURE.Simulation_large",
      "ADMIXTURE.Empirical"
    )
    
    series.labels <- c(
      "Truth (small)",
      "Truth (large)",
      "ADMIXTURE (small)",
      "ADMIXTURE (large)",
      "ADMIXTURE (empirical)"
    )
    
    series.colors <- c(
      "tspop.Simulation_small" = "#00AEDB",
      "tspop.Simulation_large" = "#007E9F",
      "ADMIXTURE.Simulation_small" = "#FA9EBC",
      "ADMIXTURE.Simulation_large" = "#E66791",
      "ADMIXTURE.Empirical" = "#B83264"
    )
    
    method.shape <- 25
    
  } else if (plot.method == "fastStructure") {
    
    series.breaks <- c(
      "tspop.Simulation_small",
      "tspop.Simulation_large",
      "fastStructure.Simulation_small",
      "fastStructure.Simulation_large",
      "fastStructure.Empirical"
    )
    
    series.labels <- c(
      "Truth (small)",
      "Truth (large)",
      "fastStructure (small)",
      "fastStructure (large)",
      "fastStructure (empirical)"
    )
    
    series.colors <- c(
      "tspop.Simulation_small" = "#00AEDB",
      "tspop.Simulation_large" = "#007E9F",
      "fastStructure.Simulation_small" = "#FFDBD1",
      "fastStructure.Simulation_large" = "#E69B88",
      "fastStructure.Empirical" = "#B9584A"
    )
    
    method.shape <- 24
    
  } else {
    stop("plot.method must be 'ADMIXTURE' or 'fastStructure'")
  }
  
  anc.length.sd.plot.df <- anc.summary %>%
    filter(
      chrom != "all",
      (
        data.type %in% c("Simulation_small", "Simulation_large") &
          method %in% c("tspop", plot.method)
      ) |
        (
          data.type == "Empirical" &
            method == plot.method
        )
    ) %>%
    group_by(
      chrom,
      pop,
      data.type,
      method
    ) %>%
    summarize(
      estimate = mean(sd, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    left_join(
      chrom.lengths,
      by = "chrom"
    ) %>%
    mutate(
      chr.len.mb = if_else(
        data.type == "Empirical",
        chr_len_after_qc,
        chr_len
      ) / 1e6,
      chrom = as.character(chrom),
      method = factor(
        method,
        levels = c("tspop", "ADMIXTURE", "fastStructure")
      ),
      data.type = factor(
        data.type,
        levels = c(
          "Simulation_small",
          "Simulation_large",
          "Empirical"
        )
      ),
      series = interaction(
        method,
        data.type,
        sep = ".",
        drop = TRUE
      ),
      label = (
        data.type == "Empirical" &
          as.integer(chrom) %% 3 == 1
      )
    )
  
  ggplot(
    anc.length.sd.plot.df,
    aes(
      x = chr.len.mb,
      y = estimate,
      fill = series,
      color = series,
      linetype = data.type,
      group = series
    )
  ) +
    geom_smooth(
      method = "lm",
      formula = y ~ x,
      se = FALSE,
      linewidth = 2
    ) +
    geom_point(
      aes(shape = method),
      color = "black",
      stroke = 1,
      size = 4
    ) +
    geom_text_repel(
      data = anc.length.sd.plot.df %>%
        filter(label),
      aes(
        x = chr.len.mb,
        y = estimate,
        label = chrom
      ),
      inherit.aes = FALSE,
      color = "black",
      size = 5,
      seed = 123,
      min.segment.length = 0,
      force = 15,
      force_pull = 5,
      box.padding = 0.5,
      point.padding = 0.25,
      max.iter = 10000,
      show.legend = FALSE
    ) +
    scale_fill_manual(
      name = NULL,
      breaks = series.breaks,
      labels = series.labels,
      values = series.colors
    ) +
    scale_color_manual(
      name = NULL,
      breaks = series.breaks,
      labels = series.labels,
      values = series.colors
    ) +
    scale_shape_manual(
      name = NULL,
      values = c(
        "tspop" = 21,
        "ADMIXTURE" = 25,
        "fastStructure" = 24
      )
    ) +
    scale_linetype_manual(
      name = NULL,
      values = c(
        "Simulation_small" = "solid",
        "Simulation_large" = "longdash",
        "Empirical" = "dashed"
      )
    ) +
    labs(
      x = "Chromosome length (Mb)",
      y = "Mean standard deviation of African ancestry",
      title = paste(
        "Chromosome Ancestry Standard Deviation versus Chromosome Length:",
        plot.method
      )
    ) +
    guides(
      color = "none",
      shape = "none",
      linetype = "none",
      fill = guide_legend(
        order = 1,
        override.aes = list(
          color = "black",
          linetype = 0,
          shape = c(
            21,
            21,
            method.shape,
            method.shape,
            method.shape
          )
        )
      )
    ) +
    theme_bw(base_size = 24) +
    theme(
      legend.position = "top",
      panel.grid.minor = element_blank()
    )
}

## plot ----
anc.length.sd.admixture.plot <- make.anc.length.sd.plot(
  plot.method = "ADMIXTURE"
)

anc.length.sd.fastStructure.plot <- make.anc.length.sd.plot(
  plot.method = "fastStructure"
)

anc.length.sd.admixture.plot
anc.length.sd.fastStructure.plot


# mean density curves ----


## calculate number of individuals per ancestry proportion bin ----
breaks <- seq(0, 1, by = 0.05)
anc.hist.chrom.levels <- c(anc.summary.chrom.levels, "all")

## Prepare ancestry values ----
anc.hist <- anc %>%
  filter(
    role == "ADX",
    chrom %in% anc.hist.chrom.levels,
    !is.na(afr.q),
    between(afr.q, 0, 1),
    method %in% c("tspop", "ADMIXTURE", "fastStructure")
  ) %>%
  transmute(
    chrom,
    rep,
    pop,
    data.type,
    method,
    ancestry = afr.q
  ) %>%
  mutate(
    chrom = factor(
      chrom,
      levels = anc.hist.chrom.levels
    ),
    data.type = factor(
      data.type,
      levels = c(
        "Simulation_small",
        "Simulation_large",
        "Empirical"
      )
    ),
    method = factor(
      method,
      levels = c(
        "tspop",
        "ADMIXTURE",
        "fastStructure"
      )
    ),
    series = interaction(
      method,
      data.type,
      sep = ".",
      drop = TRUE
    )
  )

## Histogram summary ----
anc.hist.summary <- anc.hist %>%
  group_by(
    chrom,
    rep,
    method,
    data.type,
    series
  ) %>%
  group_modify(
    ~ {
      h <- hist(
        .x$ancestry,
        breaks = breaks,
        plot = FALSE
      )
      
      tibble(
        xmin = head(h$breaks, -1),
        xmax = tail(h$breaks, -1),
        xmid = h$mids,
        count = h$counts,
        frac = h$counts / sum(h$counts)
      )
    }
  ) %>%
  ungroup() %>%
  group_by(
    chrom,
    method,
    data.type,
    series,
    xmin,
    xmax,
    xmid
  ) %>%
  summarize(
    mean.frac = mean(frac, na.rm = TRUE),
    sd.frac = if (n() > 1) {
      sd(frac, na.rm = TRUE)
    } else {
      NA_real_
    },
    n.rep = n_distinct(rep),
    .groups = "drop"
  ) %>%
  mutate(
    ymin = if_else(
      data.type == "Empirical",
      NA_real_,
      pmax(0, mean.frac - 2 * sd.frac)
    ),
    ymax = if_else(
      data.type == "Empirical",
      NA_real_,
      mean.frac + 2 * sd.frac
    )
  )

## histogram style plot ----
make.anc.hist.plot <- function(
    plot.method,
    chromosome.levels = anc.hist.chrom.levels
) {
  
  if (plot.method == "ADMIXTURE") {
    
    series.breaks <- c(
      "tspop.Simulation_small",
      "tspop.Simulation_large",
      "ADMIXTURE.Simulation_small",
      "ADMIXTURE.Simulation_large",
      "ADMIXTURE.Empirical"
    )
    
    series.labels <- c(
      "Truth (small)",
      "Truth (large)",
      "ADMIXTURE (small)",
      "ADMIXTURE (large)",
      "ADMIXTURE (empirical)"
    )
    
    series.colors <- c(
      "tspop.Simulation_small" = "#00AEDB",
      "tspop.Simulation_large" = "#007E9F",
      "ADMIXTURE.Simulation_small" = "#FA9EBC",
      "ADMIXTURE.Simulation_large" = "#E66791",
      "ADMIXTURE.Empirical" = "#B83264"
    )
    
  } else if (plot.method == "fastStructure") {
    
    series.breaks <- c(
      "tspop.Simulation_small",
      "tspop.Simulation_large",
      "fastStructure.Simulation_small",
      "fastStructure.Simulation_large",
      "fastStructure.Empirical"
    )
    
    series.labels <- c(
      "Truth (small)",
      "Truth (large)",
      "fastStructure (small)",
      "fastStructure (large)",
      "fastStructure (empirical)"
    )
    
    series.colors <- c(
      "tspop.Simulation_small" = "#00AEDB",
      "tspop.Simulation_large" = "#007E9F",
      "fastStructure.Simulation_small" = "#FFDBD1",
      "fastStructure.Simulation_large" = "#E69B88",
      "fastStructure.Empirical" = "#B9584A"
    )
    
  } else {
    stop("plot.method must be 'ADMIXTURE' or 'fastStructure'")
  }
  
  plot.df <- anc.hist.summary %>%
    filter(
      (
        method == "tspop" &
          data.type %in% c(
            "Simulation_small",
            "Simulation_large"
          )
      ) |
        method == plot.method,
      !(chrom == "all" & method == "tspop")
    )
  
  simulation.error.df <- plot.df %>%
    filter(
      data.type %in% c(
        "Simulation_small",
        "Simulation_large"
      ),
      !is.na(sd.frac)
    )
  
  dodge.position <- position_dodge(
    width = diff(breaks)[1] * 0.95
  )
  
  ggplot(
    plot.df,
    aes(
      x = xmid,
      y = mean.frac,
      fill = series,
      group = series
    )
  ) +
    geom_col(
      position = dodge.position,
      width = diff(breaks)[1] * 0.95,
      color = "black",
      linewidth = 0.5
    ) +
    geom_errorbar(
      data = simulation.error.df,
      aes(
        ymin = ymin,
        ymax = ymax
      ),
      position = dodge.position,
      width = 0.01,
      linewidth = 0.5,
      color = "black"
    ) +
    facet_wrap(
      ~ chrom,
      ncol = 3,
      drop = TRUE
    ) +
    scale_x_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.2)
    ) +
    scale_fill_manual(
      name = NULL,
      breaks = series.breaks,
      labels = series.labels,
      values = series.colors
    ) +
    labs(
      x = "African ancestry",
      y = "Mean fraction of individuals per bin",
      title = paste(
        "African Ancestry Distributions:",
        plot.method
      ),
      subtitle = paste(
        "Simulation error bars represent ±2 SD",
        "across replicate simulations"
      )
    ) +
    theme_bw(base_size = 24) +
    theme(
      legend.position = "top",
      legend.box = "horizontal",
      legend.title = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      panel.spacing = unit(1, "lines")
    )
}

anc.hist.admixture.plot <- make.anc.hist.plot(
  plot.method = "ADMIXTURE"
)

anc.hist.fastStructure.plot <- make.anc.hist.plot(
  plot.method = "fastStructure"
)

anc.hist.admixture.plot
anc.hist.fastStructure.plot

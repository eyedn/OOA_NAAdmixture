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


chrom.levels <- c(as.character(22:1), "all")
chrom.lengths <- read_table(
  "~/proj/1000GenomeNYGC_hg38_karatas/ONEKG_chr_lens.tsv", col_names = TRUE, 
  show_col_types = FALSE
  ) %>%
  mutate(chr = str_replace(as.character(chr), "chr", "")) %>%
  rename(chrom = chr)

sim.data.dir <- "~/scratch/OOA_NAAdmixture_small_cp/stats"
emp.data.dir <- "~/scratch/OOA_NAAdmixture_1kG/stats"


# read in data ----


## read in simulated data ----
sim.anc.chr <- map_dfr(
  1:22,
  \(chr) {
    file.path(sim.data.dir, glue("ancestry.chr{chr}.parquet")) %>%
      read_parquet() %>%
      mutate(
        rep = as.numeric(rep), afr_tspop = as.numeric(afr_tspop),
        eur_tspop = as.numeric(eur_tspop), afr_q = as.numeric(afr_q),
        eur_q = as.numeric(eur_q), component_1_q = NA, component_2_q = NA,
        component_3_q = NA, component_4_q = NA, component_5_q = NA,
        span = as.numeric(span), k = 0, method = "sim"
      )
    }
  )

sim.anc.genome <- read_parquet(
  file.path(sim.data.dir, glue("ancestry.parquet"))
  ) %>%
  mutate(
    rep = as.numeric(rep), chrom = "all", afr_tspop = as.numeric(afr_tspop),
    eur_tspop = as.numeric(eur_tspop), afr_q = as.numeric(afr_q),
    eur_q = as.numeric(eur_q), span = as.numeric(span), component_1_q = NA, 
    component_2_q = NA, component_3_q = NA, component_4_q = NA, 
    component_5_q = NA, k = 0, method = "sim"
    )


## read in ADMIXTURE data ----
emp.ADMIXTURE.anc.chr <- map_dfr(
  1:22,
  \(chr) {
    file.path(
      emp.data.dir, glue("ancestry_ADMIXTURE_multik.chr{chr}.parquet")
      ) %>%
      read_parquet() %>%
      mutate(
        chrom = as.character(chrom), sample_id = as.character(sample_id),
        afr_tspop = NA, eur_tspop = NA, afr_q = NA, eur_q = NA,
        method = "ADMIXTURE"
      )
    }
  )

emp.ADMIXTURE.anc.genome <- read_parquet(
  file.path(emp.data.dir, glue("ancestry_ADMIXTURE_multik.parquet"))
  ) %>%
  mutate(
    chrom = "all", sample_id = as.character(sample_id), afr_tspop = NA,
    eur_tspop = NA, afr_q = NA, eur_q = NA, method = "ADMIXTURE"
    )


## read in fastStructure data ----
emp.fastStructure.anc.chr <- map_dfr(
  1:22,
  \(chr) {
    file.path(
      emp.data.dir, glue("ancestry_fastStructure_multik.chr{chr}.parquet")
      ) %>%
      read_parquet() %>%
      mutate(
        chrom = as.character(chrom), sample_id = as.character(sample_id),
        afr_tspop = NA, eur_tspop = NA, afr_q = NA, eur_q = NA,
        method = "fastStructure"
      )
    }
  )

emp.fastStructure.anc.genome <- read_parquet(
  file.path(
    emp.data.dir, glue("ancestry_fastStructure_multik.parquet")
    )
  ) %>%
  mutate(
    chrom = "all",
    sample_id = as.character(sample_id),
    afr_tspop = NA,
    eur_tspop = NA,
    afr_q = NA,
    eur_q = NA,
    method = "fastStructure"
    )


## read in FLARE data ----
# NOTE: run analuze_anc.Rmd from traceAdmix first to get anc.lengths
emp.FLARE.anc.chr <- anc.lengths %>%
  rename(chrom = chr, vcf_sample_id = sam, component_1_q = anc1.frac) %>%
  filter(vcf_sample_id %in% emp.ADMIXTURE.anc.chr$vcf_sample_id) %>%
  mutate(
    chrom = as.character(chrom), rep = 0, sample_id = NA, k = 0,
    component_2_q = 1 - component_1_q, component_3_q = NA, component_4_q = NA,
    component_5_q = NA, span = NA, afr_tspop = NA, eur_tspop = NA,
    afr_q = NA, eur_q = NA, method = "FLARE"
    ) %>%
  select(colnames(emp.ADMIXTURE.anc.chr))


## combine ancestry data ----
anc <- rbind(
  sim.anc.chr, emp.ADMIXTURE.anc.chr, emp.fastStructure.anc.chr, 
  emp.FLARE.anc.chr, sim.anc.genome, emp.ADMIXTURE.anc.genome, 
  emp.fastStructure.anc.genome
  ) %>%
  mutate(
    role = ifelse(
      pop %in% c("AFR", "YRI"), "AFR",
      ifelse(pop %in% c("EUR", "CEU"), "EUR", "ADX")
      )
    )


## read in fastStructure chooseK results ----
emp.fastStructure.choseK.chr <- map_dfr(
  1:22,
  \(chr) {
    file.path(
      emp.data.dir, glue("fastStructure_chooseK.chr{chr}.parquet")
    ) %>%
      read_parquet()
    }
  )

emp.fastStructure.choseK.genom <- read_parquet(
  file.path(emp.data.dir, glue("fastStructure_chooseK.parquet"))
  )


## admixutre plot of optimal k for fastStructure and ADMIXTURE ----
k.use <- 2
admixPlot.df <- anc %>%
  filter(k == k.use, role == "ADX") %>%
  mutate(
    component_1_q = as.numeric(component_1_q),
    component_2_q = as.numeric(component_2_q),
    component_3_q = as.numeric(component_3_q),
    component_4_q = as.numeric(component_4_q),
    component_5_q = as.numeric(component_5_q),
    chrom = factor(chrom, levels = chrom.levels)
    ) %>%
  pivot_longer(
    starts_with("component_"), names_to = "component", values_to = "q"
    ) %>%
  filter(component %in% paste0("component_", seq_len(k.use), "_q"))

fastStructure.admixPlot <- ggplot(
  admixPlot.df %>% filter(method == "fastStructure"),
  aes(x = vcf_sample_id, y = q, fill = component)
  ) +
  geom_col() +
  facet_wrap(~chrom, ncol = 4) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  labs(
    x = NULL, y = "Ancestry proportion",
    title = "fastStructure Admixture plots"
    ) +
  theme_bw(base_size = 24) +
  theme(
    axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    strip.background = element_blank(), panel.spacing.x = unit(0.1, "lines")
    )
fastStructure.admixPlot

ADMIXTURE.admixPlot <- ggplot(
  admixPlot.df %>% filter(method == "ADMIXTURE"),
  aes(x = vcf_sample_id, y = q, fill = component)
  ) +
  geom_col() +
  facet_wrap(~chrom, ncol = 4) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  labs(
    x = NULL, y = "Ancestry proportion",
    title = "ADMIXTURE Admixture plots"
  ) +
  theme_bw(base_size = 24) +
  theme(
    axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    strip.background = element_blank(), panel.spacing.x = unit(0.1, "lines")
  )
ADMIXTURE.admixPlot


# summarize data ----
sim.anc.col <- "afr_q"

# target.component determined fom above waterfall plots
emp.anc.col.lookup <- bind_rows(
  tibble(
    method = "fastStructure",
    chrom = chrom.levels,
    target.component = c(
      "component_1_q", "component_1_q", "component_1_q", "component_2_q",
      "component_1_q", "component_2_q", "component_1_q", "component_1_q",
      "component_2_q", "component_1_q", "component_1_q", "component_1_q",
      "component_2_q", "component_2_q", "component_2_q", "component_2_q",
      "component_2_q", "component_1_q", "component_1_q", "component_2_q",
      "component_1_q", "component_2_q", "component_2_q"
    )
  ),
  tibble(
    method = "ADMIXTURE",
    chrom = chrom.levels,
    target.component = c(
      "component_2_q", "component_2_q", "component_2_q", "component_2_q",
      "component_2_q", "component_2_q", "component_2_q", "component_2_q",
      "component_2_q", "component_2_q", "component_2_q", "component_2_q",
      "component_2_q", "component_2_q", "component_2_q", "component_2_q",
      "component_2_q", "component_2_q", "component_2_q", "component_2_q",
      "component_2_q", "component_2_q", "component_2_q"
    )
  ),
  tibble(
    method = "FLARE",
    chrom = chrom.levels[1:22],
    target.component = c(
      "component_1_q", "component_1_q", "component_1_q", "component_1_q",
      "component_1_q", "component_1_q", "component_1_q", "component_1_q",
      "component_1_q", "component_1_q", "component_1_q", "component_1_q",
      "component_1_q", "component_1_q", "component_1_q", "component_1_q",
      "component_1_q", "component_1_q", "component_1_q", "component_1_q",
      "component_1_q", "component_1_q"
    )
  )
)

sim.summary <- anc %>% 
  filter(role == "ADX", method == "sim", k == 0) %>%
  group_by(chrom, rep, pop, method) %>%
  summarize(
    mean = mean(.data[[sim.anc.col]], na.rm = TRUE), 
    sd = sd(.data[[sim.anc.col]], na.rm = TRUE), 
    .groups = "drop"
    )

emp.summary <- anc %>% 
  filter(role == "ADX", method != "sim", k %in% c(0, k.use)) %>%
  pivot_longer(
    starts_with("component_"), names_to = "component", values_to = "component_q"
    ) %>%
  left_join(emp.anc.col.lookup, by = c("method", "chrom")) %>%
  filter(component == target.component) %>%
  group_by(chrom, rep, pop, method) %>%
  summarize(mean = mean(component_q), sd = sd(component_q), .groups = "drop")

anc.summary <- bind_rows(sim.summary, emp.summary)


## boxplots of mean ancestry ----
anc.summary.mean.medians <- anc.summary %>%
  group_by(chrom, pop, method) %>%
  summarize(median = median(mean), .groups = "drop") %>%
  mutate(
    y = match(chrom, chrom.levels),
    point_type = case_when(
      pop == "ADX" ~ "Simulation",
      method == "fastStructure" ~ "fastStructure",
      method == "ADMIXTURE" ~ "ADMIXTURE",
      method == "FLARE" ~ "FLARE"
    )
  )

meanAfrAnc.plot <- ggplot(
  anc.summary %>% filter(pop == "ADX"),
  aes(x = factor(chrom, levels = chrom.levels), y = mean)
  ) +
  geom_boxplot(
    aes(fill = "Simulation"), width = 0.55, linewidth = 1, outlier.shape = NA,
    color = "black",
    ) +
  geom_jitter(
    aes(fill = "Simulation"), shape = 21, color = "black", stroke = 1,
    size = 1.5, width = 0.15,height = 0
    ) +
  geom_point(
    data = anc.summary.mean.medians %>% filter(method != "sim"),
    aes(
      x = factor(chrom, levels = chrom.levels), y = median, fill = point_type,
      shape = point_type
      ),
    inherit.aes = FALSE, color = "black", stroke = 1, size = 3.5
    ) +
  scale_fill_manual( 
    name = NULL,
    breaks = c(
      "Simulation", "ADMIXTURE", "fastStructure", "FLARE"
      ),
    values = c(
      Simulation = "#00aedb", ADMIXTURE = "#FA9EBC", fastStructure = "#FFDBD1",
      FLARE = "#F8F3EA"
      )
    ) +
  scale_shape_manual(
    name = NULL,
    breaks = c("Simulation", "ADMIXTURE", "fastStructure", "FLARE"),
    values = c(Simulation = 21, ADMIXTURE = 25, fastStructure = 24, FLARE = 23)
    ) +
  labs(
    x = "Chromosome", y = "Mean African ancestry",
    title = "Mean African Ancestry of Admixed Chromosomes"
  ) +
  guides(
    fill = guide_legend(
      override.aes = list(
        shape = c(21, 25, 24, 23), color = "black", alpha = 1
        ),
      order = 1
      ),
    shape = "none"
  ) +
  theme_bw(base_size = 24) +
  theme(legend.position = "top", panel.grid.minor = element_blank())
meanAfrAnc.plot


## boxplots of sd ancestry ----
anc.summary.sd.medians <- anc.summary %>%
  group_by(chrom, pop, method) %>%
  summarize(median = median(sd), .groups = "drop") %>%
  mutate(
    y = match(chrom, chrom.levels),
    point_type = case_when(
      pop == "ADX" ~ "Simulation",
      method == "fastStructure" ~ "fastStructure",
      method == "ADMIXTURE" ~ "ADMIXTURE",
      method == "FLARE" ~ "FLARE"
    )
  )

sdAfrAnc.plot <- ggplot(
  anc.summary %>% filter(pop == "ADX"),
  aes(x = factor(chrom, levels = chrom.levels), y = sd)
  ) +
  geom_boxplot(
    aes(fill = "Simulation"), width = 0.55, linewidth = 1, outlier.shape = NA,
    color = "black"
    ) +
  geom_jitter(
    aes(fill = "Simulation"), shape = 21, color = "black", stroke = 1, 
    size = 1.5, width = 0.15, height = 0
    ) +
  geom_point(
    data = anc.summary.sd.medians %>% filter(method != "sim"),
    aes(
      x = factor(chrom, levels = chrom.levels), y = median, fill = point_type,
      shape = point_type
      ),
    inherit.aes = FALSE, color = "black", stroke = 1, size = 3.5
    ) +
  scale_fill_manual(
    name = NULL,
    breaks = c("Simulation", "ADMIXTURE", "fastStructure", "FLARE"),
    values = c(
      Simulation = "#00aedb", ADMIXTURE = "#FA9EBC", fastStructure = "#FFDBD1",
      FLARE = "#F8F3EA"
      )
    ) +
  scale_shape_manual(
    name = NULL,
    breaks = c("Simulation", "ADMIXTURE", "fastStructure", "FLARE"),
    values = c(Simulation = 21, ADMIXTURE = 25, fastStructure = 24, FLARE = 23)
    ) +
  labs(
    x = "Chromosome",
    y = "Standard deviation of African ancestry",
    title = "Standard Deviation of African Ancestry of Admixed Chromosomes",
  ) +
  guides(
    fill = guide_legend(
      override.aes = list(
        shape = c(21, 25, 24, 23), color = "black", alpha = 1
        ),
      order = 1
    ),
    shape = "none"
  ) +
  theme_bw(base_size = 24) +
  theme(legend.position = "top", panel.grid.minor = element_blank())
sdAfrAnc.plot


# mean ancestry vs chromosome length ----
anc.length.mean.summary <- anc.summary %>%
  group_by(chrom, pop, method) %>%
  summarize(estimate = median(mean), .groups = "drop") %>%
  left_join(chrom.lengths, by = "chrom") %>%
  mutate(
    chr.len.plot = ifelse(pop == "ADX", chr_len, chr_len_after_qc),
    chr.len.mb = chr.len.plot / 1e6,
    method.plot = ifelse(pop == "ADX", "Simulation", method),
    label = case_when(
      method.plot == "Simulation" ~ TRUE,
      method.plot == "ADMIXTURE" & as.integer(chrom) %% 3 == 1 ~ TRUE,
      method.plot == "fastStructure" & as.integer(chrom) %% 3 == 2 ~ TRUE,
      method.plot == "FLARE" & as.integer(chrom) %% 3 == 0 ~ TRUE,
      TRUE ~ FALSE
    )
  )


## stats ----


### sim stats ----
anc.length.mean.sim.fit <- lm(
  estimate ~ chr.len.mb,
  data = anc.length.mean.summary %>% filter(method.plot == "Simulation")
  )

summary(anc.length.mean.sim.fit)
confint(anc.length.mean.sim.fit)
# par(mfrow = c(2,2)); plot(anc.length.mean.sim.fit)


### empirical stats ----
anc.leng.emp.mean.summary <- anc.length.mean.summary %>% 
  filter(method.plot != "Simulation")

anc.length.mean.emp.fit <- lm(
  estimate ~ chr.len.mb + method.plot, data = anc.leng.emp.mean.summary
  )

summary(anc.length.mean.emp.fit)
confint(anc.length.mean.emp.fit)
# par(mfrow = c(2,2)); plot(anc.length.mean.emp.fit)


### methods stats ----
anc.length.mean.empMethods.fit <- lm(
  estimate ~ chr.len.mb * method.plot, data = anc.leng.emp.mean.summary
  )

emp.method.fits <- by(
  anc.leng.emp.mean.summary, anc.leng.emp.mean.summary$method.plot,
  function(df) lm(estimate ~ chr.len.mb, data = df)
  )

anc.length.mean.empMethods.compare <-
  anova(anc.length.mean.emp.fit, anc.length.mean.empMethods.fit)


### main stats tables and results ----
extract.lm.stats <- function(fit, dataset, method, term = "chr.len.mb") {
  td <- broom::tidy(fit)
  gl <- broom::glance(fit)
  ci <- confint(fit)

  tibble(
    dataset = dataset, method = method,
    estimate = td$estimate[td$term == "chr.len.mb"],
    conf.low = ci["chr.len.mb", 1], conf.high = ci["chr.len.mb", 2],
    p.value = td$p.value[td$term == "chr.len.mb"],
    r.squared = gl$r.squared, adj.r.squared = gl$adj.r.squared
  )
}

anc.length.mean.results <- bind_rows(
  extract.lm.stats(
    anc.length.mean.sim.fit,
    dataset = "Simulation",
    method = "Simulation"
    ),
  extract.lm.stats(
    anc.length.mean.emp.fit,
    dataset = "Empirical",
    method = "All methods"
    ),
  bind_rows(
    lapply(
      names(emp.method.fits),
      function(m) extract.lm.stats(
        emp.method.fits[[m]], dataset = "Empirical", method = m
        )
      )
    )
  )
anc.length.mean.results
anc.length.mean.empMethods.compare


## plot ----
anc.length.mean.plot <- ggplot(
  anc.length.mean.summary,
  aes(
    x = chr.len.plot / 1e6, y = estimate, fill = method.plot,
    shape = method.plot
    )
  ) +
  geom_smooth(color = "black", method = "lm", se = FALSE, linewidth = 3) +
  geom_smooth(
    aes(color = method.plot), method = "lm", se = FALSE, linewidth = 2
    ) +
  geom_point(color = "black", stroke = 1, size = 4) +
  geom_text_repel(
    data = anc.length.mean.summary %>% filter(label),
    aes(x = chr.len.plot / 1e6, y = estimate, label = chrom),
    color = "black", size = 5, seed = 123, min.segment.length = 0,
    force = 15, force_pull = 5, box.padding = 0.5, point.padding = 0.25,
    max.iter = 10000, show.legend = FALSE
    ) +
  scale_fill_manual(
    name = NULL,
    breaks = c("Simulation", "ADMIXTURE", "fastStructure", "FLARE"),
    values = c(
      Simulation = "#00aedb", ADMIXTURE = "#FA9EBC", fastStructure = "#FFDBD1",
      FLARE = "#F8F3EA"
      )
    ) +
  scale_color_manual(
    name = NULL,
    breaks = c("Simulation", "ADMIXTURE", "fastStructure", "FLARE"),
    values = c(
      Simulation = "#00aedb", ADMIXTURE = "#FA9EBC", fastStructure = "#FFDBD1",
      FLARE = "#F8F3EA"
      )
    ) +
  scale_shape_manual(
    name = NULL,
    breaks = c("Simulation", "ADMIXTURE", "fastStructure", "FLARE"),
    values = c(Simulation = 21, ADMIXTURE = 25, fastStructure = 24, FLARE = 23)
    ) +
  labs(
    x = "Chromosome length (Mb)", y = "Expected African ancestry",
    title = "Chromosome ancestry versus chromosome length"
    ) +
  guides(
    color = "none",
    shape = "none",
    fill = guide_legend(
      override.aes = list(shape = c(21, 25, 24, 23), color = "black")
      )
    ) +
  theme_bw(base_size = 24) +
  theme(legend.position = "top")
anc.length.mean.plot


# sd ancestry vs chromosome length ----
anc.length.sd.summary <- anc.summary %>%
  group_by(chrom, pop, method) %>%
  summarize(estimate = mean(sd), .groups = "drop") %>%
  left_join(chrom.lengths, by = "chrom") %>%
  mutate(
    chr.len.plot = ifelse(pop == "ADX", chr_len, chr_len_after_qc),
    chr.len.mb = chr.len.plot / 1e6,
    method.plot = ifelse(pop == "ADX", "Simulation", method),
    label = case_when(
      method.plot == "Simulation" ~ TRUE,
      method.plot == "ADMIXTURE" & as.integer(chrom) %% 3 == 1 ~ TRUE,
      method.plot == "fastStructure" & as.integer(chrom) %% 3 == 2 ~ TRUE,
      method.plot == "FLARE" & as.integer(chrom) %% 3 == 0 ~ TRUE,
      TRUE ~ FALSE
    )
  )


## stats ----


### sim stats ----
anc.length.sd.sim.fit <- lm(
  estimate ~ chr.len.mb,
  data = anc.length.sd.summary %>%
    filter(method.plot == "Simulation")
)

summary(anc.length.sd.sim.fit)
confint(anc.length.sd.sim.fit)
# par(mfrow = c(2, 2)); plot(anc.length.sd.sim.fit)


### empirical stats ----
anc.length.emp.sd.summary <- anc.length.sd.summary %>%
  filter(method.plot != "Simulation")

anc.length.sd.emp.fit <- lm(
  estimate ~ chr.len.mb + method.plot,
  data = anc.length.emp.sd.summary
)

summary(anc.length.sd.emp.fit)
confint(anc.length.sd.emp.fit)
# par(mfrow = c(2, 2)); plot(anc.length.sd.emp.fit)


### methods stats ----
anc.length.sd.empMethods.fit <- lm(
  estimate ~ chr.len.mb * method.plot,
  data = anc.length.emp.sd.summary
)

sd.method.fits <- by(
  anc.length.emp.sd.summary,
  anc.length.emp.sd.summary$method.plot,
  function(df) lm(estimate ~ chr.len.mb, data = df)
)

anc.length.sd.empMethods.compare <-
  anova(anc.length.sd.emp.fit,anc.length.sd.empMethods.fit)


### main stats tables and results ----
anc.length.sd.results <- bind_rows(
  extract.lm.stats(
    anc.length.sd.sim.fit,
    dataset = "Simulation",
    method = "Simulation"
  ),
  extract.lm.stats(
    anc.length.sd.emp.fit,
    dataset = "Empirical",
    method = "All methods"
  ),
  bind_rows(
    lapply(
      names(sd.method.fits),
      function(m)
        extract.lm.stats(
          sd.method.fits[[m]],
          dataset = "Empirical",
          method = m
        )
    )
  )
)
anc.length.sd.results
anc.length.sd.empMethods.compare

## plot ----
anc.length.sd.plot <- ggplot(
  anc.length.sd.summary,
  aes(x = chr.len.mb, y = estimate, fill = method.plot, shape = method.plot)
  ) +
  geom_smooth(color = "black", method = "lm", se = FALSE, linewidth = 3) +
  geom_smooth(
    aes(color = method.plot), method = "lm", se = FALSE, linewidth = 2
    ) +
  geom_point(color = "black", stroke = 1, size = 4) +
  geom_text_repel(
    data = anc.length.sd.summary %>% filter(label),
    aes(x = chr.len.mb, y = estimate, label = chrom), color = "black",
    size = 5, seed = 123, min.segment.length = 0, force = 15, force_pull = 5,
    box.padding = 0.5, point.padding = 0.25, max.iter = 10000, 
    show.legend = FALSE
    ) +
  scale_fill_manual(
    name = NULL,
    breaks = c("Simulation", "ADMIXTURE", "fastStructure", "FLARE"),
    values = c(
      Simulation = "#00aedb", ADMIXTURE = "#FA9EBC",fastStructure = "#FFDBD1",
      FLARE = "#F8F3EA"
      )
    ) +
  scale_color_manual(
    name = NULL,
    breaks = c("Simulation", "ADMIXTURE", "fastStructure", "FLARE"),
    values = c(
      Simulation = "#00aedb", ADMIXTURE = "#FA9EBC", fastStructure = "#FFDBD1",
      FLARE = "#F8F3EA"
      )
    ) +
  scale_shape_manual(
    name = NULL,
    breaks = c("Simulation", "ADMIXTURE", "fastStructure", "FLARE"),
    values = c(Simulation = 21, ADMIXTURE = 25, fastStructure = 24, FLARE = 23)
    ) +
  labs(
    x = "Chromosome length (Mb)", 
    y = "Expected standard deviation of African ancestry",
    title = "Chromosome ancestry standard deviation versus chromosome length"
  ) +
  guides(
    color = "none", shape = "none",
    fill = guide_legend(
      override.aes = list(shape = c(21, 25, 24, 23), color = "black")
      )
    ) +
  theme_bw(base_size = 24) +
  theme(legend.position = "top")
anc.length.sd.plot


# mean density curves ----


## calculate number of individuals per ancestry proportion bin ----
breaks <- seq(0, 1, by = 0.05)

anc.hist <- bind_rows(
  anc %>% 
    filter(
      role == "ADX", method == "sim", k == 0, !is.na(.data[[sim.anc.col]])
      ) %>%
    mutate(ancestry = .data[[sim.anc.col]]), 
  anc %>%
    filter(role == "ADX", method != "sim", k %in% c(0, 2)) %>%
    pivot_longer(
      starts_with("component_"), names_to = "component", 
      values_to = "component_q"
      ) %>%
    left_join(emp.anc.col.lookup, by = c("method", "chrom")) %>%
    filter(component == target.component) %>%
    mutate(ancestry = component_q)
  ) %>%
  mutate(
    source = case_when(
      method == "sim" ~ "Simulation",
      method == "fastStructure" ~ "fastStructure",
      method == "ADMIXTURE" ~ "ADMIXTURE",
      method == "FLARE" ~ "FLARE"
    )
  )

anc.hist.summary <- anc.hist %>%
  group_by(chrom, rep, source) %>%
  group_modify(~{
    h <- hist(.x$ancestry, breaks = breaks, plot = FALSE)
    tibble(
      xmin  = head(h$breaks, -1), xmax  = tail(h$breaks, -1), xmid  = h$mids,
      count = h$counts, frac  = h$counts / sum(h$counts)
      )
    }) %>%
  ungroup() %>%
  group_by(source, chrom, xmin, xmax, xmid) %>%
  summarize(
    mean.frac = mean(frac), sd.frac = if (n() > 1) sd(frac) else 0, 
    .groups = "drop"
    ) %>%
  mutate(
    source = factor(
      source, levels = c("Simulation", "ADMIXTURE", "fastStructure", "FLARE"
      )
    )
  ) %>%
  group_by(chrom, source) %>%
  mutate(
    ymin = pmax(0, mean.frac - 2 * sd.frac),
    ymax = mean.frac + 2 * sd.frac,
    ymin.smooth = predict(loess(ymin ~ xmid, span = 0.5)),
    ymax.smooth = predict(loess(ymax ~ xmid, span = 0.5))
  ) %>%
  ungroup()


## histogram style plot ----
anc.hist.plot <- ggplot(
  anc.hist.summary %>% 
    filter(
      # source != "Simulation",
      !(chrom == "all" & source == "Simulation")
      # chrom %in% c("1", "15", "22", "all")
      ),
  aes(x = xmid, y = mean.frac, fill = source)
  ) +
  geom_col(
    position = "dodge", width = diff(breaks)[1] * 0.95, color = "black",
    linewidth = 0.5
    ) +
  geom_errorbar(
    aes(ymin = pmax(0, mean.frac - 2 * sd.frac), ymax = mean.frac + 2 * sd.frac,
      group = source
      ),
    position = position_dodge(width = diff(breaks)[1] * 0.95),
    width = 0.025, linewidth = 0.5, color = "black"
    ) +
  facet_wrap(~factor(chrom, levels = rev(chrom.levels))) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_fill_manual(
    name = NULL, 
    breaks = c("Simulation", "ADMIXTURE", "fastStructure", "FLARE"),
    values = c(
      Simulation = "#00aedb", ADMIXTURE = "#FA9EBC", fastStructure = "#FFDBD1",
      FLARE = "#F8F3EA"
      )
    ) +
  labs(
    x = "African ancestry", y = "Mean fraction of individuals per bin",
    title = "African Ancestry Distributions of Admixed chromosomes",
    subtitle = "Error bars represent ±2 SD"
  ) +
  theme_bw(base_size = 24) +
  theme(
    legend.position = "top", legend.box = "horizontal", 
    legend.title = element_blank(), panel.grid.minor = element_blank(),
    strip.background = element_blank(), 
    strip.text = element_text(face = "bold"), panel.spacing = unit(1, "lines")
    )
anc.hist.plot


## smooth loess style plot ----
anc.smooth.plot <- ggplot(
  anc.hist.summary %>% 
    filter(
      !(chrom == "all" & source == "Simulation"),
      chrom %in% c("1", "15", "22", "all")
    ),
  aes(x = xmid, y = mean.frac)
  ) +
  geom_ribbon(
    data = anc.hist.summary %>% filter(
      source == "Simulation", 
      !(chrom == "all" & source == "Simulation"),
      chrom %in% c("1", "15", "22", "all")
      ),
    aes(ymin = ymin.smooth,ymax = ymax.smooth, fill = source), color = "black",
    alpha = 0.2
    ) +
  geom_smooth(
    aes(group = source), color = "black",method = "loess", se = FALSE, 
    span = 0.5, linewidth = 3
    ) +
  geom_smooth(
    aes(color = source), method = "loess", se = FALSE, span = 0.5, 
    linewidth = 2
    ) +
  facet_wrap(~factor(chrom, levels = rev(chrom.levels)), ncol = 2) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_color_manual(
    name = NULL, 
    breaks = c("Simulation", "ADMIXTURE", "fastStructure", "FLARE"),
    values = c(
      Simulation = "#00aedb", ADMIXTURE = "#FA9EBC", fastStructure = "#FFDBD1",
      FLARE = "#F8F3EA"
      )
    ) +
  scale_fill_manual(
    name = NULL, 
    breaks = c("Simulation", "ADMIXTURE", "fastStructure", "FLARE"),
    values = c(
      Simulation = "#00aedb", ADMIXTURE = "#FA9EBC", fastStructure = "#FFDBD1",
      FLARE = "#F8F3EA"
    )
  ) +
  labs(
    x = "African ancestry", y = "Mean fraction of individuals",
    title = "African Ancestry Distributions of Admixed chromosomes",
    subtitle = "Ribbon represent ±2 SD"
    ) +
  guides(
    fill = "none"
  ) +
  theme_bw(base_size = 24) +
  theme(
    legend.position = "top", legend.box = "horizontal", 
    legend.title = element_blank(), panel.grid.minor = element_blank(),
    strip.background = element_blank(), 
    strip.text = element_text(face = "bold"), panel.spacing = unit(1, "lines")
  )
anc.smooth.plot

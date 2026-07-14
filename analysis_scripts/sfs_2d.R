# ______________________________________________________________________________
# Aydin Karatas
# ___
# University of Southern California
# Department of Quantitative and Computational Biology 
# Mooney Lab
# ___
# sfs_2d.R
# ______________________________________________________________________________


library(tidyverse)
library(nanoparquet)
library(glue)
library(ggridges)
library(scales)
library(RColorBrewer)

# read in data ----
data.dir <- "~/scratch/OOA_NAAdmixture_small/stats"
sfs2D.chr <- map_dfr(
  1:1,
  \(chr) {
    file.path(data.dir, glue("sfs_2d.chr{chr}.parquet")) |>
      read_parquet() 
  }
)

sfs2D.genWide <- read_parquet(
  file.path(data.dir, "sfs_2d.parquet")
) %>%
  mutate(
    chrom = "all"
  )

sfs2D <- bind_rows(sfs2D.chr, sfs2D.genWide) %>%
  filter(pop1_count != "0", pop2_count != "0")

# summarize sfs ----
sfs2D.summary <- sfs2D %>%
  mutate(
    pop1_count = as.numeric(pop1_count),
    pop2_count = as.numeric(pop2_count),
    count = as.numeric(count),
    frac = count / sum(count)
  ) %>%
  group_by(chrom, pop1, pop2, pop1_count, pop2_count) %>%
  summarize(
    mean.frac = mean(frac),
    sd.frac = sd(frac),
    mean.count = mean(count),
    sd.count = sd(count),
    .groups = "drop"
  )

# plot sfs ----
max.der.allele <- 25
ggplot(
  sfs2D.summary %>%
    filter(
      chrom == "all",
      pop1_count <= max.der.allele,
      pop2_count <= max.der.allele
    ),
  aes(
    x = pop1_count,
    y = pop2_count,
    fill = mean.count
  )
) +
  geom_tile(color = "white", linewidth = 0.1) +
  geom_abline(slope = 1, intercept = 0, colour = "black", linetype = "dashed") +
  facet_wrap(
    ~(paste(pop1, pop2, sep = " vs. "))
  ) +
  scale_fill_gradientn(
    colours = rev(rainbow(10)),
    trans = "log10",
    name = "Mean\ncount"
  ) +
  coord_equal() +
  scale_x_continuous(
    breaks = seq(0,max.der.allele,5)
  ) +
  scale_y_continuous(
    breaks = seq(0,max.der.allele,5)
  ) +
  labs(
    x = "Derived allele count (Population 1)",
    y = "Derived allele count (Population 2)",
    title = "Two-dimensional Site Frequency Spectrum",
    subtitle = "Mean counts across replicates"
  ) +
  theme_bw(base_size = 24) +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  )
  
# 04_03_Cation_visualization.R
#
# Purpose:
#   Visualize habitat differences in exchangeable soil cations.
#
# Input:
#   Output/04_Soil_analysis/01_Cation_processing/
#     - All_ex_cations_long_cor_blank.rds
#     - All_soil_analysis_wilcox_test_result.csv
#
# Output:
#   Output/04_Soil_analysis/03_Cation_visualization/
#     - All_cation_plot_cor.pdf
#
# Plot:
#   - Boxplots with individual observations
#   - Faceted by site
#   - Significant habitat differences are indicated by asterisks.
#
# R version:
#   R 4.5.3
#
# Packages:
#   tidyverse
#   here
#   ggpubr

# ======================================================================
# 1. Setup
# ======================================================================
library(tidyverse)
library(here)
library(ggpubr)

# Input and output directories ------------------------------------------
input <- here("Output", "04_Soil_analysis", "01_Cation_processing")
output <- here("Output", "04_Soil_analysis", "03_Cation_visualization")
dir.create(output, showWarnings = FALSE, recursive = TRUE)

# Plot order ------------------------------------------------------------
habitat_order <- c("Solfatara field", "Forest edge")
cation_order <- c("Ca", "Mg", "K", "Na", "Al", "Mn", "Zn", "Pb", "Ni")

# ======================================================================
# 2. Load processed data
# ======================================================================
# Blank-corrected and LOD-filtered cation data.
cation_data_long <- readRDS(file.path(input, "All_ex_cations_long_cor_blank.rds")) |>
  mutate(
    habitat = factor(habitat, levels = habitat_order),
    Measurement = factor(Measurement, levels = cation_order)
  )

# Wilcoxon test results and effect sizes.
test_df <- read.csv(file.path(input, "All_soil_analysis_wilcox_test_result.csv")) |>
  mutate(Measurement = factor(Measurement, levels = cation_order))

# ======================================================================
# 3. Prepare significance annotations
# ======================================================================
sig_df <- cation_data_long |>
  group_by(site, Measurement) |>
  summarise(y_pos = max(Value, na.rm = TRUE) * 1.1, .groups = "drop") |>
  left_join(test_df |> select(site, Measurement, p.signif), by = c("site", "Measurement")) |>
  filter(!is.na(p.signif))

sig_df2 <- sig_df |>
  mutate(
    group1 = "Solfatara field",
    group2 = "Forest edge",
    xmin = as.numeric(Measurement) - 0.1875,
    xmax = as.numeric(Measurement) + 0.1875
  )

# ======================================================================
# 4. Create boxplot
# ======================================================================
p <- ggplot(cation_data_long, aes(x = Measurement, y = Value, fill = habitat)) +
  geom_boxplot(
    width = 0.65,
    outlier.shape = NA,
    color = "black",
    linewidth = 0.3,
    position = position_dodge(width = 0.75)
  ) +
  geom_jitter(
    aes(color = habitat),
    position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.75),
    size = 1.6,
    alpha = 0.6,
    show.legend = FALSE
  ) +
  # significance brackets + asterisks
  stat_pvalue_manual(
    data = sig_df2,
    label = "p.signif",
    y.position = "y_pos",
    xmin = "xmin",
    xmax = "xmax",
    size = 5,
    tip.length = 0.01,
    bracket.size = 0.4,
    inherit.aes = FALSE
  ) +
  scale_x_discrete(
    labels = c(
      "Ca" = expression(Ca^"2+"),
      "Mg" = expression(Mg^"2+"),
      "K"  = expression(K^"+"),
      "Na" = expression(Na^"+"),
      "Al" = expression(Al^"3+"),
      "Mn" = expression(Mn^"2+"),
      "Zn" = expression(Zn^"2+"),
      "Pb" = expression(Pb^"2+"),
      "Ni" = expression(Ni^"2+")
    )
  ) +
  facet_wrap(~ site, scales = "free_y", nrow = 1) +
  scale_fill_manual(
    values = c(
      "Solfatara field" = "#D55E00",
      "Forest edge" = "#0072B2"
    )
  ) +
  scale_color_grey(start = 0.35, end = 0.65) +
  labs(
    x = NULL,
    y = "The amount of exchangeable cations (mg kg-1 dry soil)",
    fill = "Habitat"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 16),
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 15),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 16),
    legend.position = "bottom"
  )

p

# ======================================================================
# 5. Save figure
# ======================================================================
ggsave(
  filename = file.path(output, "All_cation_plot_cor.pdf"),
  plot = p,
  width = 12,
  height = 6,
  device = "pdf"
)
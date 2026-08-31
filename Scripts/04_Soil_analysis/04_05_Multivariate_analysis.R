# 04_05_Multivariate_analysis.R
#
# Purpose:
#   Perform multivariate analysis (PCA, PERMANOVA, and dispersion test) 
#   on soil chemical properties (exchangeable cations and pH) to evaluate
#   environmental differences between habitats.
#
# Input:
#   Output/04_Soil_analysis/01_Cation_processing/
#     - All_ex_cations_corrected_by_blank.csv
#   Data/Soil_analysis/
#     - metadata_patch_level.csv
#
# Output:
#   Output/04_Soil_analysis/05_Multivariate_analysis/
#     - PCA_coordinates_each_sample.csv
#     - PCA_loadings.csv
#     - Cation_PCA_plot_top6_arrow.pdf
#     - PERMANOVA_results.csv
#     - Dispersion_test.csv
#
# Analysis & Plot:
#   - PCA on scaled soil cations and pH.
#   - PCA biplot with top 6 loading vectors represented as ion-labeled arrows.
#   - Stratified PERMANOVA (by site) to test habitat effects.
#   - Multivariate dispersion test (betadisper) to compare environmental variation.
#
# R version:
#   R 4.5.3
#
# Packages:
#   tidyverse
#   here
#   vegan
#   effectsize
#   ggrepel
#   broom
#   rstatix

# ======================================================================
# 1. Setup
# ======================================================================
library(here)
library(tidyverse)
library(vegan)
library(effectsize)
library(ggrepel)
library(broom)
library(rstatix)

# Input and output directories ------------------------------------------
input <- here("Output", "04_Soil_analysis", "01_Cation_processing")
output <- here("Output", "04_Soil_analysis", "05_Multivariate_analysis")
dir.create(output, showWarnings = FALSE, recursive = TRUE)

# ======================================================================
# 2. Load and merge data
# ======================================================================
exchangeable_cations <- read.csv(here(input, "All_ex_cations_corrected_by_blank.csv"))

# Remove rows with NA values (B31-Al) and the column Co (all zero values)
excatf <- na.omit(exchangeable_cations) |> select(-Co)

# Load pH metadata
meta_patch <- read.csv(here("Data", "Soil_analysis", "metadata_patch_level.csv")) |>
  select(Sample_ID, site, habitat, pH)

# Merge cations and pH data by common Sample_IDs
commons <- intersect(excatf$Sample_ID, meta_patch$Sample_ID)
excatff <- excatf[excatf$Sample_ID %in% commons, ]
pHf <- meta_patch[meta_patch$Sample_ID %in% commons, ]

merged_df <- excatff |> left_join(pHf, by = "Sample_ID")

# ======================================================================
# 3. PCA analysis
# ======================================================================
# Extract numerical environmental variables
env_vars <- merged_df |> select(-Sample_ID, -site, -habitat)

# PCA with scaled variables
pca_res <- prcomp(env_vars, scale. = TRUE)

# PCA coordinates for each sample 
pca_df <- as.data.frame(pca_res$x) |>
  mutate(
    habitat = factor(merged_df$habitat, levels = c("Solfatara field", "Forest edge")),
    site = merged_df$site
  )

write.csv(pca_df, file.path(output, "PCA_coordinates_each_sample.csv"), row.names = FALSE)

# Contributions of PC1 and PC2
explained_var <- summary(pca_res)$importance[2, 1:2] * 100
explained_var

# PCA loadings of each variable
loadings <- as.data.frame(pca_res$rotation) |>
  rownames_to_column(var = "Factors")

write.csv(loadings, file.path(output, "PCA_loadings.csv"), row.names = FALSE)

# ======================================================================
# 4. Create PCA plot
# ======================================================================
loadings <- as.data.frame(pca_res$rotation[, 1:2])
loadings$var <- rownames(loadings)

# Strength of each loading vector in PC1-PC2 space
loadings$strength <- sqrt(loadings$PC1^2 + loadings$PC2^2)

# Select top 6 variables based on loading strength
top_vars <- loadings |>
  arrange(desc(strength)) |>
  slice(1:6) |>
  pull(var)

ion_labels <- c(
  "Mn" = "Mn^{\"2+\"}",
  "Ca" = "Ca^{\"2+\"}",
  "Pb" = "Pb^{\"2+\"}",
  "K"  = "K^{\"+\"}",
  "Al" = "Al^{\"3+\"}",
  "Mg" = "Mg^{\"2+\"}"
)

loadings_sel <- loadings |>
  filter(var %in% top_vars) |>
  mutate(var_ion = dplyr::recode(var, !!!ion_labels))

# Scale loading vectors to fit the PCA score space
arrow_scale <- min(
  max(abs(pca_df$PC1)) / max(abs(loadings_sel$PC1)),
  max(abs(pca_df$PC2)) / max(abs(loadings_sel$PC2))
) * 0.7

p <- ggplot(pca_df, aes(PC1, PC2)) +
  geom_point(aes(color = habitat, shape = site), size = 5) +
  stat_ellipse(aes(group = habitat, color = habitat), linetype = 2) +
  scale_color_manual(values = c("Solfatara field" = "#D55E00", "Forest edge" = "#0072B2")) +
  geom_segment(
    data = loadings_sel,
    aes(x = 0, y = 0, xend = PC1 * arrow_scale, yend = PC2 * arrow_scale),
    arrow = arrow(length = unit(0.3, "cm")),
    linewidth = 1, color = "grey20", inherit.aes = FALSE
  ) +
  geom_text_repel(
    data = loadings_sel,
    aes(x = PC1 * arrow_scale, y = PC2 * arrow_scale, label = var_ion),
    parse = TRUE, size = 5.5, box.padding = 0.6, point.padding = 0.25,
    segment.color = "grey50", max.overlaps = Inf, inherit.aes = FALSE
  ) +
  labs(
    x = paste0("PC1 (", round(explained_var[1], 1), "%)"),
    y = paste0("PC2 (", round(explained_var[2], 1), "%)"),
    shape = "Site", color = "Habitat"
  ) +
  coord_equal() +
  theme_bw() +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14),
    legend.position = "bottom",
    legend.box = "vertical",
    panel.grid = element_blank(),
    legend.margin = margin(t = 0, r = 0, b = 0, l = 0),
    legend.box.margin = margin(t = 0, r = 0, b = 0, l = 0),
    legend.spacing.y = unit(0.1, "cm"),
    legend.key.size = unit(0.8, "cm")
  )

p

ggsave(
  filename = file.path(output, "Cation_PCA_plot_top6_arrow.pdf"),
  plot = p, width = 5.25, height = 4.5, dpi = 300, device = "pdf"
)

# ======================================================================
# 5. PERMANOVA
# ======================================================================
set.seed(1234)
perm <- 9999
env_scaled <- scale(env_vars)

adonis_res <- adonis2(
  env_scaled ~ habitat,
  data = merged_df,
  method = "euclidean",
  permutations = perm,
  strata = merged_df$site
)

permanova_results <- adonis_res |>
  as.data.frame() |>
  tibble::rownames_to_column("Factors")

print(permanova_results)

write.csv(permanova_results, file.path(output, "PERMANOVA_results.csv"), row.names = FALSE)

# ======================================================================
# 6. Dispersion test
# ======================================================================
set.seed(1234)
bd <- betadisper(dist(env_scaled), merged_df$habitat)

disp_res <- permutest(bd, permutations = perm, strata = merged_df$site)

disp_results <- disp_res$tab |>
  as.data.frame() |>
  tibble::rownames_to_column("Factors")

print(disp_results)

write.csv(disp_results, file.path(output, "Dispersion_test.csv"), row.names = FALSE)
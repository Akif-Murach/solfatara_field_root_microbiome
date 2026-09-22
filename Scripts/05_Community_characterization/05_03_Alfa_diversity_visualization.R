# 05_03_Alpha_diversity_plot.R
#
# Purpose:
#   Visualize alpha diversity (OTU richness) across sites, habitats, or sample 
#   types using boxplots with individual data points and statistical annotations.
#
# Input:
#   Output/05_Community_characterization/02_Alfa_diversity/<data_type>/<sample_type>/
#     - <data_type>_<sample_type>_richness_alpha_df.rds
#     - <data_type>_<sample_type>_richness_stat_result.csv
#
# Output:
#   Output/05_Community_characterization/02_Alfa_diversity/<data_type>/<sample_type>/Figure/
#     - <data_type>_<sample_type>_Richness.pdf
#
# Analysis & Plot:
#   - Generate boxplots overlayed with jittered points for OTU richness.
#   - Add p-value significance annotations for Root or Soil comparisons using ggpubr.
#   - Apply nested faceting (ggh4x) for Root&Soil multi-site/habitat visualization.
#
# R version:
#   R 4.5.3
#
# Packages:
#   here
#   ggplot2
#   ggpubr
#   ggh4x

# ======================================================================
# 1. Setup
# ======================================================================
library(here)
library(ggplot2)
library(ggpubr)
library(ggh4x)

# Analysis settings -----------------------------------------------------
data_type <- "Fungi"   # "Prokaryote" or "Fungi"
sample_type <- "Root"       # "Root" or "Root&Soil"

# Input and output directories ------------------------------------------
input <- here("Output", "05_Community_characterization", "02_Alfa_diversity", data_type, sample_type)

output <- here("Output", "05_Community_characterization", "02_Alfa_diversity", data_type, sample_type, "Figure")
dir.create(output, showWarnings = FALSE, recursive = TRUE)

# ======================================================================
# 2. Load data
# ======================================================================
alpha_df <- readRDS(file.path(input, paste0(data_type, "_", sample_type, "_richness_alpha_df.rds")))
stats_df <- read.csv(file.path(input, paste0(data_type, "_", sample_type, "_richness_stat_result.csv")))

# ======================================================================
# 3. Prepare plot parameters
# ======================================================================
comp_levels <- unique(as.character(alpha_df$sample_type))

data_label <- if (data_type == "Prokaryote") "Prokaryotic" else "Fungal"
y_label <- paste0(data_label, " OTU richness")

# ======================================================================
# 4. Generate boxplot
# ======================================================================

# ----------------------------------------------------------------------
# Root or Soil: Habitat comparison within each site
# ----------------------------------------------------------------------
if (length(comp_levels) == 1 && comp_levels %in% c("Root", "Soil")) {
  
  p <- ggplot(alpha_df, aes(x = habitat, y = Richness, fill = habitat)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_jitter(width = 0.10, size = 2, color = "grey50", alpha = 0.9) +
    facet_wrap(~site, scales = "free_y") +
    scale_fill_manual(
      values = c(
        "Solfatara field" = "#D55E00",
        "Forest edge"     = "#0072B2"
      )
    ) +
    stat_pvalue_manual(
      stats_df,
      label = "p.signif",
      xmin = "group1",
      xmax = "group2",
      y.position = "y.position",
      tip.length = 0.02,
      size = 6
    ) +
    labs(x = "Habitat", y = y_label)
  
  # ----------------------------------------------------------------------
  # Root&Soil: Root vs Soil comparison across sites/habitats
  # ----------------------------------------------------------------------
} else if (all(c("Root", "Soil") %in% comp_levels)) {
  
  # GLMM was performed across all samples.
  # Therefore, no statistical annotation is added to individual site x habitat facets.
  p <- ggplot(alpha_df, aes(x = sample_type, y = Richness, fill = sample_type)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_jitter(width = 0.10, size = 2, color = "grey50", alpha = 0.9) +
    ggh4x::facet_nested(~ site + habitat, scales = "free_y") +
    scale_fill_manual(
      values = c(
        "Root" = "#E6C800",
        "Soil" = "darkblue"
      )
    ) +
    labs(x = "Sample type", y = y_label)
  
} else {
  stop("Unexpected sample_type structure: ", paste(comp_levels, collapse = ", "))
}

# Apply common theme settings
p <- p +
  theme_bw(base_size = 15) +
  theme(
    axis.title.x = element_blank(),
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 16),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 13),
    legend.position = "bottom"
  )

# Display plot in console
print(p)

# ======================================================================
# 5. Save figure
# ======================================================================
ggsave(
  filename = file.path(output, paste0(data_type, "_", sample_type, "_Richness.pdf")),
  plot = p,
  width = 6,
  height = 6,
  units = "in",
  device = "pdf",
  dpi = 300
)
# 06_04_PCoA.R
#
# Purpose:
#   Perform Principal Coordinates Analysis (PCoA) based on Sørensen dissimilarity
#   and visualize community composition differences with 95% confidence ellipses.
#
# Input:
#   Loaded via "Scripts/06_Community_analysis/06_00_Setup.R":
#     - dist: Distance matrix (Sørensen / Bray-Curtis dissimilarity)
#     - metadata: Sample metadata
#     - data_type: "Prokaryote" or "Fungi"
#     - sample_type: "Root", "Soil", or "Root&Soil"
#
# Output:
#   Output/06_Community_analysis/04_PCoA/<data_type>/
#     - <data_type>_Root_PCoA_sorensen_th3.pdf (If Root)
#     - <data_type>_Root&Soil_PCoA_sorensen_th3.pdf (If Root&Soil or Soil)
#
# Analysis & Plot:
#   - Ordination analysis using ape::pcoa.
#   - Plotting PCoA1 and PCoA2 scores with 95% confidence ellipses (stat_ellipse).
#   - Custom shapes using ggstar::geom_star (Root) or ggplot2::geom_point (Root&Soil).
#
# R version:
#   R 4.5.3
#
# Packages:
#   here
#   ape
#   ggplot2
#   ggstar
#   dplyr
#   tibble

# ======================================================================
# 1. Setup
# ======================================================================
source(here("Scripts", "06_Community_analysis", "06_00_Setup.R"))

# Output directory ------------------------------------------------------
output <- here("Output", "06_Community_analysis", "04_PCoA", data_type)
dir.create(output, showWarnings = FALSE, recursive = TRUE)

# ======================================================================
# 2. PCoA calculation
# ======================================================================
pcoa_res <- ape::pcoa(dist)

pcoa_score <- as.data.frame(pcoa_res$vectors[, 1:2]) |>
  rownames_to_column(var = "Sample_ID")

colnames(pcoa_score) <- c("Sample_ID", "PCoA1", "PCoA2")

# Calculate variance explained (%) -------------------------------------
explain <- pcoa_res$values$Relative_eig * 100
x_label <- paste0("PCoA1 (", round(explain[1], 1), "%)")
y_label <- paste0("PCoA2 (", round(explain[2], 1), "%)")

# ======================================================================
# 3. PCoA plot
# ======================================================================

# ----------------------------------------------------------------------
# Root samples: Habitat and host plant visualization
# ----------------------------------------------------------------------
if (sample_type == "Root") {
  
  metadata$habitat <- factor(metadata$habitat)
  metadata$host <- factor(metadata$host)
  
  pcoa_df <- metadata |>
    left_join(pcoa_score, by = "Sample_ID")
  
  # Factor levels and aesthetics ---------------------------------------
  habitat_order <- c("Solfatara field", "Forest edge")
  pcoa_df$habitat <- factor(pcoa_df$habitat, levels = habitat_order)
  
  host_order <- c(
    "Enkianthus campanulatus",
    "Eubotryoides grayana",
    "Rhododendron spp.",
    "Vaccinium smallii",
    "Betula ermanii",
    "Pinus parviflora"
  )
  pcoa_df$host <- factor(pcoa_df$host, levels = host_order)
  
  star_shapes <- c(
    "Enkianthus campanulatus" = 13,
    "Eubotryoides grayana"    = 15,
    "Rhododendron spp."       = 11,
    "Vaccinium smallii"       = 12,
    "Betula ermanii"          = 5,
    "Pinus parviflora"        = 14
  )
  
  ellipse_colors <- c(
    "Solfatara field" = "#D55E00",
    "Forest edge"     = "#0072B2"
  )
  
  # Generate plot ------------------------------------------------------
  pcoa <- ggplot(
    pcoa_df,
    aes(
      x = PCoA1,
      y = PCoA2,
      fill = habitat,
      color = habitat,
      group = habitat
    )
  ) +
    geom_star(
      aes(starshape = host),
      size = 6,
      color = "black"
    ) +
    stat_ellipse(
      aes(color = habitat),
      type = "norm",
      linewidth = 0.7,
      linetype = "dashed"
    ) +
    scale_fill_manual(values = ellipse_colors) +
    scale_color_manual(values = ellipse_colors, guide = "none") +
    scale_starshape_manual(values = star_shapes) +
    theme_bw(base_size = 15) +
    theme(
      panel.grid = element_blank(),
      panel.background = element_rect(fill = "white"),
      axis.text = element_text(size = 15),
      legend.title = element_text(size = 15)
    ) +
    guides(
      fill = guide_legend(
        title = "Habitat",
        override.aes = list(starshape = 13, alpha = 1)
      ),
      color = "none",
      starshape = guide_legend(
        title = "Host plant identity",
        override.aes = list(fill = "black")
      )
    ) +
    labs(x = x_label, y = y_label)
  
  # Save plot ----------------------------------------------------------
  print(pcoa)
  
  ggsave(
    filename = here(output, paste0(data_type, "_Root_PCoA_sorensen_th3.pdf")),
    plot = pcoa,
    width = 10,
    height = 7,
    units = "in",
    device = "pdf",
    dpi = 300
  )
  
  # ----------------------------------------------------------------------
  # Root + Soil samples: Habitat and sample_type visualization
  # ----------------------------------------------------------------------
} else {
  
  metadata$sample_type <- factor(metadata$sample_type)
  metadata$habitat <- factor(metadata$habitat)
  
  pcoa_df <- metadata |>
    left_join(pcoa_score, by = "Sample_ID")
  
  # Factor levels and aesthetics ---------------------------------------
  habitat_order <- c("Solfatara field", "Forest edge")
  pcoa_df$habitat <- factor(pcoa_df$habitat, levels = habitat_order)
  
  sample_type_order <- c("Root", "Soil")
  pcoa_df$sample_type <- factor(pcoa_df$sample_type, levels = sample_type_order)
  
  ellipse_colors <- c(
    "Root" = "#E6C800",
    "Soil" = "darkblue"
  )
  
  shape_values <- c(21, 24)
  
  # Generate plot ------------------------------------------------------
  pcoa <- ggplot(pcoa_df, aes(x = PCoA1, y = PCoA2)) +
    geom_point(
      aes(fill = sample_type, shape = habitat),
      size = 6,
      stroke = 0.7,
      color = "black"
    ) +
    stat_ellipse(
      aes(color = sample_type, group = sample_type),
      type = "norm",
      linewidth = 0.7,
      linetype = "dashed"
    ) +
    scale_fill_manual(values = ellipse_colors) +
    scale_color_manual(values = ellipse_colors, guide = "none") +
    scale_shape_manual(values = shape_values) +
    theme_bw(base_size = 15) +
    theme(
      panel.grid = element_blank(),
      panel.background = element_rect(fill = "white"),
      axis.text = element_text(size = 15),
      legend.title = element_text(size = 15)
    ) +
    guides(
      fill = guide_legend(
        title = "Sample type",
        override.aes = list(shape = 21)
      ),
      shape = guide_legend(
        title = "Habitat",
        override.aes = list(fill = "black")
      )
    ) +
    labs(x = x_label, y = y_label)
  
  # Save plot ----------------------------------------------------------
  print(pcoa)
  
  ggsave(
    filename = here(output, paste0(data_type, "_Root&Soil_PCoA_sorensen_th3.pdf")),
    plot = pcoa,
    width = 9,
    height = 7,
    units = "in",
    device = "pdf",
    dpi = 300
  )
}
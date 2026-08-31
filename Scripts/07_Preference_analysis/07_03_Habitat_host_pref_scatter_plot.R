# 07_03_Habitat_host_pref_scatter_plot.R
#
# Purpose:
#   Visualize the relationship between habitat preference and host preference 
#   for fungal or prokaryotic OTUs via scatter plot.
#
# Input:
#   Loaded via "Scripts/07_Preference_analysis/07_00_Setup.R":
#     - data_type: Target data group ("fungi", "prokaryote", etc.)
#     - threshold: Filtering threshold value
#   Files read:
#     - Output/07_Preference_analysis/2DP/<data_type>/habitat/th<threshold>/2DP_Zvalue_habitat.rds
#     - Output/07_Preference_analysis/2DP/<data_type>/habitat/th<threshold>/2DP_two_sided_FDR_habitat.rds
#     - Output/07_Preference_analysis/dprime/<data_type>/microbe/th<threshold>/dprime_Zvalue_microbe.rds
#     - Output/07_Preference_analysis/dprime/<data_type>/microbe/th<threshold>/dprime_two_sided_FDR_microbe.rds
#     - Data/<data_type>/Seqdata/OTU_merge_taxonomylist.rds
#
# Output:
#   - Output/07_Preference_analysis/Proc_data/<data_type>/th<threshold>/preference_data.rds
#   - Output/07_Preference_analysis/Figures/<data_type>/th<threshold>/<data_type>_habitat_host_pref_scatter.pdf
#
# Analysis & Visualization:
#   1. Load habitat- and host-preference results.
#   2. Retain OTUs with complete preference information.
#   3. Classify OTUs according to significant habitat and host preferences.
#   4. Add taxonomic annotations (highest available rank + OTU_ID).
#   5. Select representative OTUs for text labeling (top host, top habitat, both significant).
#   6. Generate and save the habitat-host preference scatter plot.
#
# R version:
#   R 4.5.3
#
# Packages:
#   here
#   tidyverse
#   ggrepel
#   ggtext

# ======================================================================
# 1. Setup & Package Load
# ======================================================================
library(here)
library(tidyverse)
library(ggrepel)
library(ggtext)

# Load analysis settings -----------------------------------------------
source(here("Scripts", "07_Preference_analysis", "07_00_Setup.R"))

# Define input/output directory paths ----------------------------------
pref_dir <- here("Output", "07_Preference_analysis")

habitat_pref <- file.path(pref_dir, "2DP", data_type, "habitat", paste0("th", threshold))
host_pref    <- file.path(pref_dir, "dprime", data_type, "microbe", paste0("th", threshold))

seq_dir <- here("Data", data_type, "Seqdata")

output  <- file.path(pref_dir, "Proc_data", data_type, paste0("th", threshold))
output2 <- file.path(pref_dir, "Figures", data_type, paste0("th", threshold))

dir.create(output, showWarnings = FALSE, recursive = TRUE)
dir.create(output2, showWarnings = FALSE, recursive = TRUE)

# ======================================================================
# 2. Load Preference Data
# ======================================================================
# Habitat preference ---------------------------------------------------
habitat_z <- readRDS(file.path(habitat_pref, "2DP_Zvalue_habitat.rds")) |>
  as.data.frame()

habitat_p <- readRDS(file.path(habitat_pref, "2DP_two_sided_FDR_habitat.rds")) |>
  as.data.frame()

# Host preference ------------------------------------------------------
host_z <- readRDS(file.path(host_pref, "dprime_Zvalue_microbe.rds")) |>
  as.data.frame()

host_p <- readRDS(file.path(host_pref, "dprime_two_sided_FDR_microbe.rds")) |>
  as.data.frame()

# ======================================================================
# 3. Merge & Save Preference Data
# ======================================================================
common_otus <- Reduce(
  intersect,
  list(
    rownames(habitat_z),
    rownames(habitat_p),
    rownames(host_z),
    rownames(host_p)
  )
)

preference_df <- tibble(
  OTU_ID             = common_otus,
  habitat_preference = habitat_z[common_otus, 1],
  habitat_p          = habitat_p[common_otus, 1],
  host_preference    = host_z[common_otus, 1],
  host_p             = host_p[common_otus, 1]
)

saveRDS(preference_df, file.path(output, "preference_data.rds"))

# ======================================================================
# 4. Classify OTUs by Significant Preference
# ======================================================================
preference_df <- preference_df |>
  mutate(
    Significant = case_when(
      host_p < 0.05 & habitat_p < 0.05 ~ "Habitat&Host",
      host_p < 0.05                    ~ "Host",
      habitat_p < 0.05                 ~ "Habitat",
      TRUE                             ~ "Non"
    ),
    Siglabel = case_when(
      Significant == "Habitat&Host" & habitat_preference > 0  ~ "Solfatara field&Host",
      Significant == "Habitat&Host" & habitat_preference <= 0 ~ "Forest edge&Host",
      Significant == "Habitat"      & habitat_preference > 0  ~ "Solfatara field",
      Significant == "Habitat"      & habitat_preference <= 0 ~ "Forest edge",
      Significant == "Host"                                   ~ "Host",
      TRUE                                                    ~ "Non"
    )
  )

# ======================================================================
# 5. Add Taxonomic Labels
# ======================================================================
tax <- readRDS(file.path(seq_dir, "OTU_merge_taxonomylist.rds")) |>
  as.data.frame() |>
  rownames_to_column("OTU_ID")

tax_label <- tax |>
  mutate(
    across(
      c(Species, Genus, Family, Order, Class, Phylum),
      ~ ifelse(
        is.na(.) | grepl("unidentified|Incertae_sedis", ., ignore.case = TRUE),
        NA_character_,
        .
      )
    ),
    label = coalesce(Genus, Family, Order, Class, Phylum, "Unidentified"),
    label = paste0(label, "\n[", OTU_ID, "]")
  ) |>
  select(OTU_ID, label)

plot_data <- preference_df |>
  left_join(tax_label, by = "OTU_ID")

# ======================================================================
# 6. Select OTUs to be Labeled in Plot
# ======================================================================
# Top host-preference OTUs ---------------------------------------------
top_host <- plot_data |>
  filter(Significant == "Host") |>
  arrange(desc(host_preference)) |>
  slice_head(n = 3)

# Top habitat-preference OTUs ------------------------------------------
if (sum(plot_data$Siglabel == "Solfatara field", na.rm = TRUE) > 0) {
  top_habitat_positive <- plot_data |>
    filter(Significant == "Habitat") |>
    arrange(desc(habitat_preference)) |>
    slice_head(n = 2)
  
  top_habitat_negative <- plot_data |>
    filter(Significant == "Habitat") |>
    arrange(habitat_preference) |>
    slice_head(n = 3)
  
  top_habitat <- bind_rows(top_habitat_positive, top_habitat_negative)
} else {
  top_habitat <- plot_data |>
    filter(Significant == "Habitat") |>
    arrange(habitat_preference) |>
    slice_head(n = 5)
}

# Both significant OTUs ------------------------------------------------
both_significant <- plot_data |>
  filter(Significant == "Habitat&Host")

# Combine unique label targets -----------------------------------------
target_labels <- c(
  top_host$label,
  top_habitat$label,
  both_significant$label
) |>
  unique()

plot_data <- plot_data |>
  mutate(
    plot_label = if_else(label %in% target_labels, label, NA_character_)
  )

# ======================================================================
# 7. Color Palette & Factor Levels Definition
# ======================================================================
okabe_ito <- c(
  "Solfatara field"      = "#D55E00",
  "Forest edge"          = "#0072B2",
  "Host"                 = "#009E73",
  "Solfatara field&Host" = "#CC79A7",
  "Forest edge&Host"     = "#56B4E9",
  "Non"                  = "darkgray"
)

preference_levels <- c(
  "Solfatara field",
  "Forest edge",
  "Host",
  "Solfatara field&Host",
  "Forest edge&Host",
  "Non"
)

# ======================================================================
# 8. Generate Scatter Plot
# ======================================================================
p <- ggplot(
  plot_data,
  aes(
    x = host_preference,
    y = habitat_preference
  )
) +
  geom_point(
    aes(fill = Siglabel),
    size = 6,
    shape = 21,
    stroke = 0.2,
    color = "white"
  ) +
  geom_text_repel(
    aes(label = plot_label),
    vjust = -1.5,
    size = 6,
    box.padding = 0.7,
    max.overlaps = 30
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dotted",
    color = "black",
    linewidth = 0.8
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dotted",
    color = "black",
    linewidth = 0.8
  ) +
  labs(
    x = "Preference for host plants (z-standardized d')",
    y = "Preference for solfatara fields (2DP)",
    fill = "Significance"
  ) +
  scale_fill_manual(
    values = okabe_ito,
    limits = preference_levels,
    labels = c(
      "Solfatara-\nfield",
      "Forest-\nedge",
      "Host",
      "Solfatara-\nfield&Host",
      "Forest-\nedge&Host",
      "Non"
    ),
    drop = FALSE
  ) +
  theme_bw() +
  theme(
    axis.title.x  = ggtext::element_markdown(size = 16),
    axis.title.y  = ggtext::element_markdown(size = 16),
    axis.text.x   = element_text(size = 16, face = "bold"),
    axis.text.y   = element_text(size = 16, face = "bold"),
    legend.title  = element_text(size = 12),
    legend.text   = element_text(size = 12),
    plot.subtitle = element_textbox_simple()
  )

# Display plot ---------------------------------------------------------
p

# ======================================================================
# 9. Save Figure
# ======================================================================
ggsave(
  filename = file.path(
    output2,
    paste0(data_type, "_habitat_host_pref_scatter.pdf")
  ),
  plot = p,
  width = 9,
  height = 7,
  units = "in",
  dpi = 300,
  device = "pdf"
)
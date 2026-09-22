# 07_05_Host_pref_heatmap.R
#
# Purpose:
#   Visualize host preference by combining:
#   1. Microbial d' preference for host plants
#   2. Host x microbe two-dimensional preference (2DP)
#   3. Host plant d' preference
#
# Input:
#   Loaded via "Scripts/07_Preference_analysis/07_00_Setup.R":
#     - data_type: Target data group ("Prokaryote" or "Fungi")
#     - threshold: Filtering threshold value
#   Files read:
#     - Host x microbe 2DP Z-scores and FDR-adjusted P-values
#     - Host d' Z-scores and FDR-adjusted P-values
#     - Microbial d' Z-scores and FDR-adjusted P-values
#     - Data/<data_type>/Seqdata/OTU_merge_taxonomylist.rds
#
# Output:
#   - Output/07_Preference_analysis/<data_type>/Figures/Host_<data_type>_pref_heatmap_th<threshold>.pdf
#
# R version:
#   R 4.5.3
#
# Packages:
#   here
#   tidyverse
#   patchwork
#   stats

# ======================================================================
# 1. Setup & Package Load
# ======================================================================
library(here)
library(tidyverse)
library(patchwork)

# Load analysis settings -----------------------------------------------
source(here("Scripts", "07_Preference_analysis", "07_00_Setup.R"))

# Define input/output directory paths ----------------------------------
pref_dir <- here("Output", "07_Preference_analysis")
seq_dir  <- here("Data", data_type, "Seqdata")

output_dir <- file.path(pref_dir, "Figures", data_type, paste0("th",threshold))
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ======================================================================
# 2. Load Preference Data
# ======================================================================
# Host x microbe two-dimensional preference (2DP) ---------------------
tdp_dir <- file.path(pref_dir, "2DP", data_type, "host", paste0("th", threshold))

tdp_z <- readRDS(file.path(tdp_dir, "2DP_Zvalue_host.rds"))
tdp_p <- readRDS(file.path(tdp_dir, "2DP_two_sided_FDR_host.rds"))

# Host d' --------------------------------------------------------------
hdp_dir <- file.path(pref_dir, "dprime", data_type, "host", paste0("th", threshold))

hdp_z <- readRDS(file.path(hdp_dir, "dprime_Zvalue_host.rds"))
hdp_p <- readRDS(file.path(hdp_dir, "dprime_two_sided_FDR_host.rds"))

# Microbial d' ---------------------------------------------------------
mdp_dir <- file.path(pref_dir, "dprime", data_type, "microbe", paste0("th", threshold))

mdp_z <- readRDS(file.path(mdp_dir, "dprime_Zvalue_microbe.rds")) |>
  as.data.frame() |>
  setNames("dprime") |>
  rownames_to_column("OTU")

mdp_p <- readRDS(file.path(mdp_dir, "dprime_two_sided_FDR_microbe.rds")) |>
  as.data.frame() |>
  setNames("p_val") |>
  rownames_to_column("OTU")

# ======================================================================
# 3. Preprocessing
# ======================================================================
# Keep OTUs even when some or all 2DP cells are untestable.
valid_otu <- rownames(tdp_z)

tdp_z <- tdp_z[valid_otu, , drop = FALSE]
tdp_p <- tdp_p[valid_otu, , drop = FALSE]
tdp_p[is.na(tdp_z)] <- NA_real_

mdp_z <- mdp_z |> filter(OTU %in% valid_otu)
mdp_p <- mdp_p |> filter(OTU %in% valid_otu)

# Select top 50 OTUs based on absolute microbial d' -------------------
top50_otu <- mdp_z |>
  arrange(desc(abs(dprime))) |>
  slice_head(n = 50) |>
  pull(OTU)

tdp_z <- tdp_z[top50_otu, , drop = FALSE]
tdp_p <- tdp_p[top50_otu, , drop = FALSE]

mdp_z <- mdp_z |> filter(OTU %in% top50_otu)
mdp_p <- mdp_p |> filter(OTU %in% top50_otu)

# ======================================================================
# 4. Taxonomy & OTU Labels
# ======================================================================
tax <- readRDS(file.path(seq_dir, "OTU_merge_taxonomylist.rds")) |>
  as.data.frame()

tax_df <- tax |>
  mutate(OTU_ID = sub("^X_", "F_", rownames(tax)))

# Create taxonomic labels for plotting ---------------------------------
taxa_table <- tax_df |>
  mutate(
    lowest_label = case_when(
      !is.na(Genus)  & !grepl("unidentified|Incertae_sedis", Genus, ignore.case = TRUE)  ~ Genus,
      !is.na(Family) & !grepl("unidentified|Incertae_sedis", Family, ignore.case = TRUE) ~ Family,
      !is.na(Order)  & !grepl("unidentified|Incertae_sedis", Order, ignore.case = TRUE)  ~ Order,
      !is.na(Class)  & !grepl("unidentified|Incertae_sedis", Class, ignore.case = TRUE)  ~ Class,
      !is.na(Phylum) & !grepl("unidentified|Incertae_sedis", Phylum, ignore.case = TRUE) ~ Phylum,
      TRUE ~ "Unidentified"
    ),
    order_label = replace_na(Order, "Unidentified"),
    label       = paste0(lowest_label, "_(", order_label, ")_[", OTU_ID, "]")
  ) |>
  select(OTU_ID, label)

# ======================================================================
# 5. Taxonomic Hierarchy & OTU Ordering
# ======================================================================
tax_hier <- tax_df |>
  select(OTU_ID, Phylum, Class, Order, Family, Genus) |>
  filter(OTU_ID %in% rownames(tdp_z)) |>
  mutate(across(c(Phylum, Class, Order, Family, Genus), ~ replace_na(.x, "Unidentified")))

# Hierarchical clustering based on 2DP Z-scores ------------------------
# Preserve the existing zero fill for ordering only; plotted Z-scores stay NA.
tdp_matrix <- tdp_z
tdp_matrix[is.na(tdp_matrix)] <- 0

otu_clustering <- hclust(dist(tdp_matrix), method = "ward.D2")

clust_df <- tibble(
  OTU_ID     = otu_clustering$labels[otu_clustering$order],
  clust_rank = seq_along(otu_clustering$order)
)

# Combine taxonomy and clustering information -------------------------
otu_order <- tax_hier |>
  left_join(clust_df, by = "OTU_ID") |>
  arrange(Phylum, Class, Order, Family, Genus, clust_rank) |>
  pull(OTU_ID)

otu_labels <- taxa_table$label[match(otu_order, taxa_table$OTU_ID)]

# Define host plant display order --------------------------------------
host_order <- c(
  "Enkianthus campanulatus",
  "Eubotryoides grayana",
  "Rhododendron spp.",
  "Vaccinium smallii",
  "Betula ermanii",
  "Pinus parviflora"
)

# ======================================================================
# 6. Build Heatmap Components
# ======================================================================
# Microbial d' heatmap -------------------------------------------------
df_tile_dprime <- mdp_z |>
  left_join(mdp_p, by = "OTU") |>
  mutate(
    OTU_label = factor(OTU, levels = otu_order, labels = otu_labels))|>
  mutate(p_val = if_else(is.na(dprime), NA_real_, p_val))

d_plt_dprime <- ggplot(
  df_tile_dprime,
  aes(x = data_type, y = OTU_label, fill = dprime)
) +
  geom_tile(width = 1, height = 1) +
  geom_text(
    aes(
      label = case_when(
        p_val < 0.001 ~ "***",
        p_val < 0.01  ~ "**",
        p_val < 0.05  ~ "*",
        TRUE          ~ ""
      )
    ),
    size = 6,
    vjust = 0.9,
    color = "white"
  ) +
  scale_fill_gradient2(
    low = "#56B4E9",
    mid = "white",
    high = "firebrick",
    midpoint = 0,
    na.value = "grey75",
    name = "Preference for host"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 10.5),
    axis.title  = element_blank()
  )

# Host x microbe 2DP heatmap -------------------------------------------
df_tile_2dp <- as.data.frame(tdp_z) |>
  rownames_to_column("OTU") |>
  pivot_longer(cols = -OTU, names_to = "Host", values_to = "Zscore") |>
  left_join(
    as.data.frame(tdp_p) |>
      rownames_to_column("OTU") |>
      pivot_longer(cols = -OTU, names_to = "Host", values_to = "p_val"),
    by = c("OTU", "Host")
  ) |>
  mutate(
    Host      = fct_relevel(Host, host_order),
    OTU_label = factor(OTU, levels = otu_order, labels = otu_labels))|>
  mutate(p_val = if_else(is.na(Zscore), NA_real_, p_val))

d_plt_2dp <- ggplot(
  df_tile_2dp,
  aes(x = Host, y = OTU_label, fill = Zscore)
) +
  geom_tile(width = 1, height = 1) +
  geom_text(
    aes(
      label = case_when(
        p_val < 0.001 ~ "***",
        p_val < 0.01  ~ "**",
        p_val < 0.05  ~ "*",
        TRUE          ~ ""
      )
    ),
    size = 6,
    vjust = 0.9,
    color = "white"
  ) +
  scale_fill_gradient2(
    low = "#0072B2",
    mid = "white",
    high = "#D55E00",
    midpoint = 0,
    na.value = "grey75",
    name = paste0("Host x ", data_type, " 2DP")
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x  = element_blank(),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title   = element_blank()
  )

# Host d' heatmap ------------------------------------------------------
df_host_dprime <- tibble(
  Host   = names(hdp_z),
  dprime = as.numeric(hdp_z),
  p_val  = as.numeric(hdp_p[names(hdp_z)])
) |>
  mutate(Host = fct_relevel(Host, host_order))|>
  mutate(p_val = if_else(is.na(dprime), NA_real_, p_val))

d_plt_host_dprime <- ggplot(
  df_host_dprime,
  aes(x = Host, y = "dprime", fill = dprime)
) +
  geom_tile(width = 1, height = 1) +
  geom_text(
    aes(
      label = case_when(
        p_val < 0.001 ~ "***",
        p_val < 0.01  ~ "**",
        p_val < 0.05  ~ "*",
        TRUE          ~ ""
      )
    ),
    size = 6,
    color = "white"
  ) +
  scale_fill_gradient2(
    low = "#56B",
    mid = "white",
    high = "darkgreen",
    midpoint = 0,
    na.value = "grey75",
    name = paste0("Preference for ", data_type)
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(size = 12, angle = 30, hjust = 1, face = "italic"),
    axis.text.y = element_blank(),
    axis.title  = element_blank()
  )

# ======================================================================
# 7. Combine & Save Heatmaps
# ======================================================================
plot_design <- "
23
#1
"

combined_plot <- d_plt_host_dprime +
  d_plt_dprime +
  d_plt_2dp +
  plot_layout(
    design  = plot_design,
    guides  = "collect",
    widths  = c(15, 85),
    heights = c(100, 5)
  )

# Display plot ---------------------------------------------------------
combined_plot

# Save figure ----------------------------------------------------------
ggsave(
  filename = file.path(
    output_dir,
    paste0("Host_", data_type, "_pref_heatmap_th", threshold, ".pdf")
  ),
  plot   = combined_plot,
  width  = ifelse(data_type == "Prokaryote", 10, 9),
  height = 12,
  units  = "in",
  dpi    = 300
)
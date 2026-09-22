# 05_01_Barplot.R
#
# Purpose:
#   Visualize microbial taxonomic composition (e.g., Order-level) across 
#   sites and habitats using barplots, and calculate relative abundances.
#
# Input:
#   Output/03_Data_filtering/Seqdata/<data_type>/<sample_type>/
#     - <data_type>_<sample_type>_coverage_rared_th1fil.rds
#   Output/03_Data_filtering/Metadata/<data_type>/<sample_type>/
#     - <data_type>_<sample_type>_metadata_th1fil.csv
#   Data/<data_type>/Seqdata/
#     - OTU_merge_taxonomylist.rds
#
# Output:
#   Output/05_Community_characterization/02_Barplot/<data_type>/<sample_type>/
#     - <data_type>_<sample_type>_<class>_prop.csv
#     - <data_type>_<sample_type>_<class>_barplot.pdf
#
# Analysis & Plot:
#   - Calculate relative abundance per taxon at the specified taxonomic rank.
#   - Order samples using Bray-Curtis distance and hierarchical clustering.
#   - Generate stacked barplots faceted by site, habitat, and sample type.
#
# R version:
#   R 4.5.3
#
# Packages:
#   tidyverse
#   ggplot2
#   vegan
#   here
#   ggh4x

# ======================================================================
# 1. Setup
# ======================================================================
library(tidyverse)
library(ggplot2)
library(vegan)
library(here)

# Load external function ------------------------------------------------
source(here("Function", "Taxa.mat.R"))

# Analysis settings -----------------------------------------------------
data_type <- "Prokaryote"      # "Prokaryote" or "Fungi"
sample_type <- "Root"     # "Root" or "Root&Soil"
class <- "Order"
limit <- 18

# Input and output directories ------------------------------------------
inputs <- here("Output", "03_Data_filtering", "Seqdata", data_type, sample_type)
inputm <- here("Output", "03_Data_filtering", "Metadata", data_type, sample_type)

output <- here("Output", "05_Community_characterization", "02_Barplot", data_type, sample_type)
dir.create(output, showWarnings = FALSE, recursive = TRUE)

# ======================================================================
# 2. Load data
# ======================================================================
# Read rarefied OTU table and metadata
seqdata <- readRDS(here(inputs, paste0(data_type, "_", sample_type, "_coverage_rared_th1fil.rds")))
metadata <- read.csv(here(inputm, paste0(data_type, "_", sample_type, "_metadata_th1fil.csv")))

tax <- readRDS(here("Data", data_type, "Seqdata", "OTU_merge_taxonomylist.rds")) |>
  as.data.frame()

# ======================================================================
# 3. Process taxonomic matrix
# ======================================================================
trr <- Taxa.mat(seqdata, tax, class)

top <- setdiff(colnames(trr)[order(colSums(trr), decreasing = TRUE)], "Unidentified")[1:limit]
print(top)

trr2 <- cbind(
  trr[, top],
  Others = rowSums(trr[, !colnames(trr) %in% c(top, "Unidentified")]),
  Unidentified = trr[, "Unidentified"]
)

trr3 <- rownames_to_column(as.data.frame(trr2), var = "sample")
trr4 <- left_join(trr3, metadata, by = join_by(sample == Sample_ID))

gtrr <- pivot_longer(
  trr4,
  cols = colnames(trr2),
  names_to = "taxa",
  values_to = "read_count"
)

# ======================================================================
# 4. Calculate relative abundance
# ======================================================================
if (sample_type != "Root&Soil") {
  taxaprop_all <- gtrr |>
    group_by(site, habitat, taxa) |>
    summarize(total_read_count = sum(read_count, na.rm = TRUE), .groups = "drop") |>
    group_by(site, habitat) |>
    mutate(relative_abundance_percent = total_read_count / sum(total_read_count) * 100) |>
    arrange(site, habitat, desc(relative_abundance_percent)) |>
    ungroup()
} else {
  taxaprop_all <- gtrr |>
    group_by(site, habitat, sample_type, taxa) |>
    summarize(total_read_count = sum(read_count, na.rm = TRUE), .groups = "drop") |>
    group_by(site, habitat, sample_type) |>
    mutate(relative_abundance_percent = total_read_count / sum(total_read_count) * 100) |>
    arrange(site, habitat, sample_type, desc(relative_abundance_percent)) |>
    ungroup()
}

# Split and list by site and habitat
genprop <- taxaprop_all |>
  group_split(site, habitat)

# Final output table
taxaprop_final <- taxaprop_all |>
  mutate(relative_abundance_percent = round(relative_abundance_percent, 2))

write.csv(
  taxaprop_final,
  file = here(output, paste0(data_type, "_", sample_type, "_", class, "_prop.csv")),
  row.names = FALSE
)

# ======================================================================
# 5. Sorting and clustering for plotting
# ======================================================================
# Determine sample order using hierarchical clustering
d <- vegdist(trr2 / rowSums(trr2), "bray")
cl <- hclust(d, method = "average")

gtrr$sample2 <- factor(gtrr$sample, levels = rownames(trr2)[cl$order])

# Color palette and factor ordering
keep_orders <- c("Others", "Unidentified", rev(top))
habitat_orders <- c("Solfatara field", "Forest edge")

if (data_type == "Prokaryote") {
  order_palette <- setNames(
    c(
      "gray50", "gray70", "#AA8BCC", "#9EC3B7", "#049CC6", "#8DA0CB",
      "#56B4E9", "#E45A1C", "#7CAE56", "#7570B3", "#FED90E", "#1B9E10",
      "#F564E3", "#A6761D", "#FF7F00", "#E31A1C", "#1F78B4", "#7F3C8D",
      "#1B9E90", "#00005B"
    ),
    keep_orders
  )
} else {
  order_palette <- setNames(
    c(
      "gray50", "gray70", "cornsilk2", "lightpink", "bisque3", "skyblue3",
      "olivedrab3", "orchid1", "orange2", "salmon3", "mediumorchid",
      "darkseagreen4", "indianred2", "dodgerblue2", "mediumpurple",
      "springgreen4", "firebrick3", "royalblue3", "firebrick4", "darkblue"
    ),
    keep_orders
  )
}

# Assign factor levels for plot formatting
gtrr$Order <- ifelse(gtrr$taxa %in% keep_orders, as.character(gtrr$taxa), "Others")
gtrr$Order <- factor(gtrr$Order, levels = keep_orders)
gtrr$habitat <- factor(gtrr$habitat, levels = habitat_orders)

# ======================================================================
# 6. Create barplot
# ======================================================================
g <- ggplot(gtrr, aes(x = sample2, y = read_count)) +
  geom_bar(aes(fill = Order, color = Order), stat = "identity", position = "fill", width = 1) +
  labs(x = "Samples", y = "Proportion of sequencing reads (%)") +
  scale_y_continuous(labels = function(x) x * 100) +
  scale_fill_manual(values = order_palette, drop = TRUE) +
  scale_color_manual(values = order_palette, drop = TRUE) +
  theme_bw(base_size = 15) +
  theme(
    axis.text.x = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(size = 16),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 13),
    legend.position = "bottom"
  )

if (sample_type != "Root&Soil") {
  g <- g + ggh4x::facet_nested(~ site + habitat, scales = "free")
} else {
  g <- g + ggh4x::facet_nested(~ site + habitat + sample_type, scales = "free")
}

print(g)

# ======================================================================
# 7. Save figure
# ======================================================================
ggsave(
  filename = here(output, paste0(data_type, "_", sample_type, "_", class, "_barplot.pdf")),
  plot = g,
  device = "pdf",
  dpi = 300,
  width = 12,
  height = 6,
  units = "in"
)
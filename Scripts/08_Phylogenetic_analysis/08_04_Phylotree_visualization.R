# 08_04_Phylotree_visualization.R
# Purpose: Visualize phylogenetic relationships of Hyaloscyphaceae OTUs with habitat & host preference.
# Workflow: 1. Load data -> 2. Format bootstrap -> 3. Match OTUs -> 4. Prepare coords -> 5. Set offsets -> 6-9. Build tree plot -> 10. Save figure
# Input: Output/05_Phylogenetic_analysis/Phylo_data/
# Output: Output/05_Phylogenetic_analysis/Fig/

#===========================================================
# Packages & Setup
#===========================================================
library(here)
library(tidyverse)
library(ape)
library(ggtree)
library(ggnewscale)
library(treeio)

input_dir <- here("Output", "08_Phylogenetic_analysis", "Phylo_data")
input_dir2 <- here("Output", "07_Preference_analysis", "Proc_data", "Fungi", "th3")
output_dir <- here("Output", "08_Phylogenetic_analysis", "Figure")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

#===========================================================
# 1. Load phylogenetic tree and trait data
#===========================================================
tree <- read.newick(file.path(input_dir, "Newick_Pairwise_0.6.nwk"))
trait_data <- readRDS(file.path(input_dir2, "preference_data.rds")) |> as.data.frame()

#===========================================================
# 2. Format bootstrap support values
#===========================================================
if (!is.null(tree$node.label)) {
  bootstrap_values <- suppressWarnings(as.numeric(tree$node.label))
  bootstrap_values <- ifelse(!is.na(bootstrap_values) & max(bootstrap_values, na.rm = TRUE) > 1, bootstrap_values / 100, bootstrap_values)
  tree$node.label <- ifelse(!is.na(bootstrap_values), round(bootstrap_values * 100), "")
}

#===========================================================
# 3. Match OTUs between tree and trait data
#===========================================================
tree_otus <- tibble(tip.label = tree$tip.label) |>
  mutate(OTU_ID = str_extract(tip.label, "F_[0-9]+"))

otu_scores <- tree_otus |> left_join(trait_data, by = "OTU_ID")

#===========================================================
# 4. Prepare tree coordinates & marker positions
#===========================================================
p_tree_base <- ggtree(tree)
tree_data <- p_tree_base$data
x_max <- max(tree_data$x, na.rm = TRUE)

tip_data <- tree_data |> filter(isTip) |> select(label, x, y)
plot_data <- otu_scores |> left_join(tip_data, by = c("tip.label" = "label"))
marker_data <- plot_data |> filter(!is.na(OTU_ID))

habitat_offset <- x_max * 0.27
host_offset <- x_max * 0.30

#===========================================================
# 5. Build Base Tree & Add Preference Markers
#===========================================================
p_tree <- ggtree(tree) +
  geom_tiplab(size = 4, offset = x_max * 0.003) +
  geom_text2(aes(label = label), 
             data = \(df) df[!df$isTip & df$label != "", , drop = FALSE],
             size = 3, hjust = 1.3, vjust = -0.5, color = "grey30") +
  theme_tree2()

# Add habitat-preference markers
p_final <- p_tree +
  geom_point(data = marker_data, aes(x = x + habitat_offset, y = y, 
                                     fill = habitat_preference), 
             shape = 22, size = 7.2, stroke = 0.15) +
  geom_text(data = marker_data, aes(x = x + habitat_offset, y = y,
                                    label = case_when(habitat_p < 0.001 ~ "***", 
                                                      habitat_p < 0.01 ~ "**", 
                                                      habitat_p < 0.05 ~ "*", TRUE ~ "")),
            size = 3.5, vjust = 1.2, color = "white") +
  scale_fill_gradient2(low = "#0072B2", mid = "white", high = "#D55E00", midpoint = 0, 
                       name = "Solfatara\nfield\npreference\n(2DP)") +
  ggnewscale::new_scale_fill()

# Add host-preference markers & Finalize
p_final <- p_final +
  geom_point(data = marker_data, aes(x = x + host_offset, y = y, fill = host_preference), 
             shape = 21, size = 7.2, stroke = 0.15) +
  geom_text(data = marker_data, aes(x = x + host_offset, y = y,
                                    label = case_when(host_p < 0.001 ~ "***", 
                                                      host_p < 0.01 ~ "**", 
                                                      host_p < 0.05 ~ "*", TRUE ~ "")),
            size = 3.5, vjust = 0.8, color = "white") +
  scale_fill_gradient2(low = "#56B4E9", mid = "white", high = "firebrick", midpoint = 0, 
                       name = "Host\npreference\n(d′)") +
  coord_cartesian(xlim = c(0, x_max * 1.25)) +
  theme(legend.position = "right", legend.box = "vertical")

# Display plot
p_final

#===========================================================
# 6. Save figure
#===========================================================
ggsave(filename = file.path(output_dir, "Phylotree_marker_pairwise0.6_th3.pdf"), 
       plot = p_final, width = 13, height = 11)
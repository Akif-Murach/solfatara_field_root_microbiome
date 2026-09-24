# 06_01_PERMANOVA.R
#
# Purpose:
#   Perform Permutational Multivariate Analysis of Variance (PERMANOVA)
#   to test community structure differences across experimental factors.
#
# Description:
#   Root samples:
#     - Test habitat x host using sequential sums of squares (by = "terms").
#     - Permutations are restricted within site.
#
#   Root & Soil samples:
#     - Test habitat effects separately within Root and Soil.
#       Permutations are restricted within site.
#     - Test Root vs Soil separately within Solfatara field and Forest edge.
#       Permutations are restricted within SID.
#     - BH FDR correction is applied separately to:
#         (1) the two habitat tests
#         (2) the two sample-type tests
#
# Input:
#   Loaded via "Scripts/06_Community_analysis/06_00_Setup.R":
#     - dist: Distance matrix
#     - metadata: Sample metadata
#     - data_type: "Prokaryote" or "Fungi"
#     - sample_type: "Root", "Soil", or "Root&Soil"
#
# Output:
#   Root:
#     Output/06_Community_analysis/01_PERMANOVA/<data_type>/
#       - <data_type>_Root_PERMANOVA_th3.csv
#
#   Root&Soil:
#     Output/06_Community_analysis/01_PERMANOVA/<data_type>/
#       - <data_type>_Root&Soil_PERMANOVA_habitat_root_th3.csv
#       - <data_type>_Root&Soil_PERMANOVA_habitat_soil_th3.csv
#       - <data_type>_Root&Soil_PERMANOVA_Solfatara_field_th3.csv
#       - <data_type>_Root&Soil_PERMANOVA_sample_type_Forest_edge_th3.csv
#
# R version: R 4.5.3
# Packages: here, vegan, dplyr, permute

# ======================================================================
# 1. Setup
# ======================================================================

library(here)
library(dplyr)
library(permute)

source(here("Scripts", "06_Community_analysis", "06_00_Setup.R"))

output <- here("Output", "06_Community_analysis", "01_PERMANOVA", data_type)
dir.create(output, showWarnings = FALSE, recursive = TRUE)

set.seed(1234)
n <- 9999

# ======================================================================
# 2. Helper functions
# ======================================================================

run_permanova <- function(dist, metadata, terms, block, n = 9999) {
  ctrl <- how(nperm = n)
  setBlocks(ctrl) <- metadata[[block]]
  formula <- reformulate(terms, response = "dist")
  vegan::adonis2(formula, data = metadata, by = "terms", permutations = ctrl) |>
    as.data.frame()
}

run_subgroup_permanova <- function(dist, metadata, filter_col, filter_val, 
                                   terms, block, n = 9999) {
  sub_meta <- metadata[metadata[[filter_col]] == filter_val, ]
  sub_dist <- as.dist(as.matrix(dist)[rownames(sub_meta), rownames(sub_meta)])
  run_permanova(dist = sub_dist, metadata = sub_meta, terms = terms, 
                block = block, n = n)
}

apply_fdr <- function(res1, res2) {
  p_vals <- c(res1$`Pr(>F)`[1], res2$`Pr(>F)`[1])
  fdr_vals <- p.adjust(p_vals, method = "BH")
  
  res1$FDR <- NA_real_; res1[1, "FDR"] <- fdr_vals[1]
  res2$FDR <- NA_real_; res2[1, "FDR"] <- fdr_vals[2]
  
  list(res1 = res1, res2 = res2)
}

# ======================================================================
# 3. Root samples: Habitat and host effects
# ======================================================================

if (sample_type == "Root") {
  permanova <- run_permanova(
    dist = dist,
    metadata = metadata,
    terms = c("habitat", "host", "habitat:host"),
    block = "site",
    n = n
  )
  
  print(permanova)
  
  out_file <- file.path(output, paste0(data_type, "_", sample_type, 
                                       "_PERMANOVA_th3.csv"))
  write.csv(permanova, file = out_file)
}

# ======================================================================
# 4. Root + Soil samples: Habitat and sample_type effects
# ======================================================================

if (sample_type == "Root&Soil") {
  
  # --- Habitat effect (permutations within site) ---
  res_habitat_root <- run_subgroup_permanova(dist, metadata, "sample_type",
                                             "Root", "habitat", "site", n)
  res_habitat_soil <- run_subgroup_permanova(dist, metadata, "sample_type", 
                                             "Soil", "habitat", "site", n)
  
  adj_habitat <- apply_fdr(res_habitat_root, res_habitat_soil)
  res_habitat_root <- adj_habitat$res1
  res_habitat_soil <- adj_habitat$res2
  
  # --- Sample type effect (permutations within SID) ---
  res_comp_sf <- run_subgroup_permanova(dist, metadata, "habitat", 
                                        "Solfatara field", "sample_type", "SID", n)
  res_comp_fe <- run_subgroup_permanova(dist, metadata, "habitat", 
                                        "Forest edge", "sample_type", "SID", n)
  
  adj_comp <- apply_fdr(res_comp_sf, res_comp_fe)
  res_comp_sf <- adj_comp$res1
  res_comp_fe <- adj_comp$res2
  
  # --- Inspect outputs ---
  print(res_habitat_root)
  print(res_habitat_soil)
  print(res_comp_sf)
  print(res_comp_fe)
  
  # --- Save results ---
  save_list <- list(
    habitat_root = res_habitat_root,
    habitat_soil = res_habitat_soil,
    sample_type_solfatara_field = res_comp_sf,
    sample_type_forest_edge = res_comp_fe
  )
  
  for (suffix in names(save_list)) {
    out_file <- file.path(output, paste0(data_type, "_", sample_type, 
                                         "_PERMANOVA_", suffix, "_th3.csv"))
    write.csv(save_list[[suffix]], file = out_file)
  }
}

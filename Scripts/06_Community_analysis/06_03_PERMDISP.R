# 06_03_PERMDISP.R
#
# Purpose:
#   Test homogeneity of multivariate dispersion across groups
#   (habitat, host, or sample_type) using permutation tests
#   (betadisper/permutest) with restricted permutations based on
#   site or sampling position (SID), depending on the comparison.
#
# Input:
#   Loaded via "Scripts/06_Community_analysis/06_00_Setup.R":
#     - dist: Distance matrix (Sorensen)
#     - metadata: Sample metadata
#     - data_type: "Prokaryote" or "Fungi"
#     - sample_type: "Root", "Soil", or "Root&Soil"
#
# Output:
#   Output/06_Community_analysis/03_PERMDISP/<data_type>/<sample_type>/
#
#   Root:
#     - <data_type>_Root_PERMDISP_habitat_th3.csv
#     - <data_type>_Root_PERMDISP_host_th3.csv
#     - <data_type>_Root_PERMDISP_host_pairwise_th3.csv
#
#   Root&Soil:
#     - <data_type>_Root&Soil_PERMDISP_habitat_root_th3.csv
#     - <data_type>_Root&Soil_PERMDISP_habitat_soil_th3.csv
#     - <data_type>_Root&Soil_PERMDISP_Solfatara_field_th3.csv
#     - <data_type>_Root&Soil_PERMDISP_sample_type_Forest_edge_th3.csv
#
# Analysis:
#   Root:
#     - Habitat: permutations restricted within site
#     - Host: permutations restricted within site
#     - Pairwise host comparisons: permutations restricted within site
#
#   Root&Soil:
#     - Habitat within Root subset: permutations restricted within site
#     - Habitat within Soil: permutations restricted within site
#     - Root vs Soil within Solfatara field: permutations restricted within SID
#     - Root vs Soil within Forest edge: permutations restricted within SID
#
#   Bias adjustment:
#     - bias.adjust = TRUE is applied to all PERMDISP analyses.
#
#   Multiple testing:
#     - BH FDR correction is applied separately to:
#         (1) the two habitat tests in the Root&Soil dataset
#         (2) the two sample-type tests in the Root&Soil dataset
#         (3) pairwise host comparisons in the Root dataset
#
# R version:
#   R 4.5.3
#
# Packages:
#   here, vegan, dplyr, permute

# ======================================================================
# 1. Setup
# ======================================================================

library(here)
library(vegan)
library(dplyr)
library(permute)

source(here("Scripts", "06_Community_analysis", "06_00_Setup.R"))

# Output directory ------------------------------------------------------
output <- here("Output", "06_Community_analysis", "03_PERMDISP", 
               data_type, sample_type)
dir.create(output, showWarnings = FALSE, recursive = TRUE)

# Random seed & permutation settings -----------------------------------
set.seed(1234)
n <- 9999

# ======================================================================
# 2. Helper function
# ======================================================================

run_betadisper <- function(dist, metadata, group, block, n = 9999, 
                           bias_adjust = TRUE) {
  # Restricted permutation design
  ctrl <- permute::how(nperm = n)
  permute::setBlocks(ctrl) <- metadata[[block]]
  
  # Multivariate dispersion & permutation test
  bd <- vegan::betadisper(dist, group, bias.adjust = bias_adjust)
  perm <- vegan::permutest(bd, permutations = ctrl)
  
  as.data.frame(perm$tab)
}

# ======================================================================
# 3. Root samples: Habitat and host effects
# ======================================================================

if (sample_type == "Root") {
  
  # Habitat effect (permutations restricted within site) ----------------
  res_habitat <- run_betadisper(
    dist = dist, metadata = metadata, group = metadata$habitat,
    block = "site", n = n, bias_adjust = TRUE
  )
  
  # Host effect (permutations restricted within site) -------------------
  res_host <- run_betadisper(
    dist = dist, metadata = metadata, group = metadata$host,
    block = "site", n = n, bias_adjust = TRUE
  )
  
  # Pairwise host comparisons (permutations restricted within site) -----
  ctrl_host <- permute::how(nperm = n)
  permute::setBlocks(ctrl_host) <- metadata$site
  
  bd_host <- vegan::betadisper(dist, metadata$host, bias.adjust = TRUE)
  perm_host_pair <- vegan::permutest(bd_host, permutations = ctrl_host, 
                                     pairwise = TRUE)
  
  pair_host <- data.frame(
    comparison = names(perm_host_pair$pairwise$observed),
    p_value = perm_host_pair$pairwise$observed,
    p_perm = perm_host_pair$pairwise$permuted
  ) |>
    mutate(FDR = p.adjust(p_perm, method = "BH"))
  
  # Save results --------------------------------------------------------
  write.csv(res_habitat, here(output, paste0(data_type, "_", sample_type,
                                             "_PERMDISP_habitat_th3.csv")))
  write.csv(res_host, here(output, paste0(data_type, "_", sample_type, 
                                          "_PERMDISP_host_th3.csv")))
  write.csv(pair_host, here(output, paste0(data_type, "_", sample_type, 
                                           "_PERMDISP_host_pairwise_th3.csv")), 
            row.names = FALSE)
}

# ======================================================================
# 4. Root + Soil samples: Habitat and sample-type effects
# ======================================================================

if (sample_type == "Root&Soil") {
  
  # Habitat effect: Root subset (restricted within site) ---------------
  info_root <- metadata |> filter(sample_type == "Root")
  dist_root <- as.dist(as.matrix(dist)[rownames(info_root), rownames(info_root)])
  
  res_habitat_root <- run_betadisper(
    dist = dist_root, metadata = info_root, group = info_root$habitat,
    block = "site", n = n, bias_adjust = TRUE
  )
  
  # Habitat effect: Soil (restricted within site) -----------------------
  info_soil <- metadata |> filter(sample_type == "Soil")
  dist_soil <- as.dist(as.matrix(dist)[rownames(info_soil), rownames(info_soil)])
  
  res_habitat_soil <- run_betadisper(
    dist = dist_soil, metadata = info_soil, group = info_soil$habitat, 
    block = "site", n = n, bias_adjust = TRUE   )    
  # FDR correction: Habitat effects ------------------------------------   
  habitat_p <- c(res_habitat_root$`Pr(>F)`[1], res_habitat_soil$`Pr(>F)`[1])
  habitat_fdr <- p.adjust(habitat_p, method = "BH")
  
  res_habitat_root$FDR <- NA_real_
  res_habitat_soil$FDR <- NA_real_
  res_habitat_root["Groups", "FDR"] <- habitat_fdr[1]
  res_habitat_soil["Groups", "FDR"] <- habitat_fdr[2]
  
  # Sample-type effect: Solfatara field (Root vs Soil, within SID) ------
  info_sf <- metadata |> filter(habitat == "Solfatara field")
  dist_sf <- as.dist(as.matrix(dist)[rownames(info_sf), rownames(info_sf)])
  
  res_comp_sf <- run_betadisper(
    dist = dist_sf, metadata = info_sf, group = info_sf$sample_type,
    block = "SID", n = n, bias_adjust = TRUE
  )
  
  # Sample-type effect: Forest edge (Root vs Soil, within SID) ---------
  info_fe <- metadata |> filter(habitat == "Forest edge")
  dist_fe <- as.dist(as.matrix(dist)[rownames(info_fe), rownames(info_fe)])
  
  res_comp_fe <- run_betadisper(
    dist = dist_fe, metadata = info_fe, group = info_fe$sample_type,    
    block = "SID", n = n, bias_adjust = TRUE   )   
  # FDR correction: Sample-type effects --------------------------------   
  comp_p <- c(res_comp_sf$`Pr(>F)`[1], res_comp_fe$`Pr(>F)`[1])
  comp_fdr <- p.adjust(comp_p, method = "BH")
  
  res_comp_sf$FDR <- NA_real_
  res_comp_fe$FDR <- NA_real_
  res_comp_sf["Groups", "FDR"] <- comp_fdr[1]
  res_comp_fe["Groups", "FDR"] <- comp_fdr[2]
  
  # Save results --------------------------------------------------------
  write.csv(res_habitat_root, here(output, paste0(data_type, "_", sample_type, 
                                                  "_PERMDISP_habitat_root_th3.csv")))
  write.csv(res_habitat_soil, here(output, paste0(data_type, "_", sample_type, 
                                                  "_PERMDISP_habitat_soil_th3.csv")))
  write.csv(res_comp_sf, here(output, paste0(data_type, "_", sample_type, 
                                             "_PERMDISP_sample_type_Solfatara_field_th3.csv")))
  write.csv(res_comp_fe, here(output, paste0(data_type, "_", sample_type, 
                                             "_PERMDISP_sample_type_Forest_edge_th3.csv")))
}
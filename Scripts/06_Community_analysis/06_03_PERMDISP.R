# 06_03_PERMDISP.R
#
# Purpose:
#   Test homogeneity of multivariate dispersion across groups (habitat, host, or sample_type)
#   using permutation tests (betadisper/permutest) with permutations constrained within sites.
#
# Input:
#   Loaded via "Scripts/06_Community_analysis/06_00_Setup.R":
#     - dist: Distance matrix (Sorensen)
#     - metadata: Sample metadata
#     - data_type: "Prokaryote" or "Fungi"
#     - sample_type: "Root", "Soil", or "Root&Soil"
#
# Output:
#   Output/03_Community_analysis/03_PERMDISP/<data_type>/<sample_type>/
#     - <data_type>_<sample_type>_PERMDISP_habitat_th3.csv (If Root)
#     - <data_type>_<sample_type>_PERMDISP_host_th3.csv (If Root)
#     - <data_type>_<sample_type>_PERMDISP_host_pairwise_th3.csv (If Root)
#     - <data_type>_<sample_type>_PERMDISP_habitat_root_th3.csv (If Root&Soil)
#     - <data_type>_<sample_type>_PERMDISP_habitat_soil_th3.csv (If Root&Soil)
#     - <data_type>_<sample_type>_PERMDISP_Solfatara_field_th3.csv (If Root&Soil)
#     - <data_type>_<sample_type>_PERMDISP_sample_type_Forest_edge_th3.csv (If Root&Soil)
#
# Analysis:
#   - Multivariate dispersion analysis via vegan::betadisper and vegan::permutest.
#   - Constrain permutations by site blocking using vegan::how.
#   - Apply BH FDR adjustment for multiple hypothesis comparisons.
#
# R version:
#   R 4.5.3
#
# Packages:
#   here
#   vegan
#   dplyr

# ======================================================================
# 1. Setup
# ======================================================================
source(here("Scripts", "06_Community_analysis", "06_00_Setup.R"))

# Output directory ------------------------------------------------------
output <- here("Output", "06_Community_analysis", "03_PERMDISP",
               data_type, sample_type)
dir.create(output, showWarnings = FALSE, recursive = TRUE)

# Random seed & permutation settings -----------------------------------
set.seed(1234)
n <- 9999

# Helper function for betadisper analysis -------------------------------
run_betadisper <- function(dist, metadata, group, n = 9999) {
  ctrl <- how(nperm = n)
  setBlocks(ctrl) <- metadata$site
  
  bd <- betadisper(dist, group)
  perm <- permutest(bd, permutations = ctrl)
  as.data.frame(perm$tab)
}

# ======================================================================
# 2. Root samples: Habitat and host effects
# ======================================================================
if (sample_type == "Root") {
  
  # Habitat effect ------------------------------------------------------
  res_habitat <- run_betadisper(
    dist = dist,
    metadata = metadata,
    group = metadata$habitat,
    n = n
  )
  
  # Host effect ---------------------------------------------------------
  res_host <- run_betadisper(
    dist = dist,
    metadata = metadata,
    group = metadata$host,
    n = n
  )
  
  # Pairwise host comparisons -------------------------------------------
  ctrl_host <- how(nperm = n)
  setBlocks(ctrl_host) <- metadata$site
  
  bd_host <- betadisper(dist, metadata$host)
  
  perm_host_pair <- permutest(
    bd_host,
    permutations = ctrl_host,
    pairwise = TRUE
  )
  
  pair_host <- data.frame(
    comparison = names(perm_host_pair$pairwise$observed),
    p_value = perm_host_pair$pairwise$observed,
    p_perm = perm_host_pair$pairwise$permuted
  ) |>
    mutate(FDR = p.adjust(p_perm, method = "BH"))
  
  # Save results --------------------------------------------------------
  write.csv(
    res_habitat,
    here(output, paste0(data_type, "_", sample_type, "_PERMDISP_habitat_th3.csv"))
  )
  
  write.csv(
    res_host,
    here(output, paste0(data_type, "_", sample_type, "_PERMDISP_host_th3.csv"))
  )
  
  write.csv(
    pair_host,
    here(output, paste0(data_type, "_", sample_type, "_PERMDISP_host_pairwise_th3.csv")),
    row.names = FALSE
  )
}

# ======================================================================
# 3. Root + Soil samples: Habitat and sample_type effects
# ======================================================================
if (sample_type == "Root&Soil") {
  
  # Habitat effect: Root ------------------------------------------------
  info_root <- metadata |>
    filter(sample_type == "Root")
  
  dist_root <- as.dist(as.matrix(dist)[rownames(info_root), rownames(info_root)])
  
  res_habitat_root <- run_betadisper(
    dist = dist_root,
    metadata = info_root,
    group = info_root$habitat,
    n = n
  )
  
  # Habitat effect: Soil ------------------------------------------------
  info_soil <- metadata |>
    filter(sample_type == "Soil")
  
  dist_soil <- as.dist(as.matrix(dist)[rownames(info_soil), rownames(info_soil)])
  
  res_habitat_soil <- run_betadisper(
    dist = dist_soil,
    metadata = info_soil,
    group = info_soil$habitat,
    n = n
  )
  
  # FDR correction for Habitat effect ----------------------------------
  habitat_p <- c(
    res_habitat_root$`Pr(>F)`[1],
    res_habitat_soil$`Pr(>F)`[1]
  )
  
  habitat_fdr <- p.adjust(habitat_p, method = "BH")
  
  res_habitat_root$FDR <- NA
  res_habitat_soil$FDR <- NA
  
  res_habitat_root["Groups", "FDR"] <- habitat_fdr[1]
  res_habitat_soil["Groups", "FDR"] <- habitat_fdr[2]
  
  # Sample type effect: Solfatara field ---------------------------------
  info_sf <- metadata |>
    filter(habitat == "Solfatara field")
  
  dist_sf <- as.dist(as.matrix(dist)[rownames(info_sf), rownames(info_sf)])
  
  res_comp_sf <- run_betadisper(
    dist = dist_sf,
    metadata = info_sf,
    group = info_sf$sample_type,
    n = n
  )
  
  # Sample type effect: Forest edge ------------------------------------
  info_fe <- metadata |>
    filter(habitat == "Forest edge")
  
  dist_fe <- as.dist(as.matrix(dist)[rownames(info_fe), rownames(info_fe)])
  
  res_comp_fe <- run_betadisper(
    dist = dist_fe,
    metadata = info_fe,
    group = info_fe$sample_type,
    n = n
  )
  
  # FDR correction for Sample type effect --------------------------------
  comp_p <- c(
    res_comp_sf$`Pr(>F)`[1],
    res_comp_fe$`Pr(>F)`[1]
  )
  
  comp_fdr <- p.adjust(comp_p, method = "BH")
  
  res_comp_sf$FDR <- NA
  res_comp_fe$FDR <- NA
  
  res_comp_sf["Groups", "FDR"] <- comp_fdr[1]
  res_comp_fe["Groups", "FDR"] <- comp_fdr[2]
  
  # Save results --------------------------------------------------------
  write.csv(
    res_habitat_root,
    here(output, paste0(data_type, "_", sample_type, "_PERMDISP_habitat_root_th3.csv"))
  )
  
  write.csv(
    res_habitat_soil,
    here(output, paste0(data_type, "_", sample_type, "_PERMDISP_habitat_soil_th3.csv"))
  )
  
  write.csv(
    res_comp_sf,
    here(output, paste0(data_type, "_", sample_type, "_PERMDISP_Solfatara_field_th3.csv"))
  )
  
  write.csv(
    res_comp_fe,
    here(output, paste0(data_type, "_", sample_type, "_PERMDISP_sample_type_Forest_edge_th3.csv"))
  )
}
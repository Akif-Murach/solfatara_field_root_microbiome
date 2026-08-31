# 06_01_PERMANOVA.R
#
# Purpose:
#   Perform Permutational Multivariate Analysis of Variance (PERMANOVA) 
#   to test community structure differences across experimental factors.
#
# Description:
#   Tests habitat x host interaction for Root samples, and habitat x sample_type 
#   interaction for Root & Soil samples. Sequential sums of squares (by = "terms") 
#   are used to account for hierarchical factor structures (habitat > host / sample_type).
#   Permutations are stratified by 'site' (strata = metadata$site).
#
# Input:
#   Loaded via "Scripts/06_Community_analysis/06_00_Setup.R":
#     - dist: Distance matrix
#     - metadata: Sample metadata
#     - data_type: "Prokaryote" or "Fungi"
#     - sample_type: "Root", "Soil", or "Root&Soil"
#
# Output:
#   Output/06_Community_analysis/01_PERMANOVA/<data_type>/
#     - <data_type>_<sample_type>_PERMANOVA_th3.csv
#
# Analysis:
#   - PERMANOVA using vegan::adonis2.
#
# R version:
#   R 4.5.3
#
# Packages:
#   here
#   vegan

# ======================================================================
# 1. Setup
# ======================================================================
source(here("Scripts", "06_Community_analysis", "06_00_Setup.R"))

# Output directory ------------------------------------------------------
output <- here("Output", "06_Community_analysis", "01_PERMANOVA", data_type)
dir.create(output, showWarnings = FALSE, recursive = TRUE)

# ======================================================================
# 2. PERMANOVA calculation
# ======================================================================
set.seed(1234)
n <- 9999

permanova <- if (sample_type == "Root") {
  vegan::adonis2(
    dist ~ habitat * host,
    data = metadata,
    by = "terms",
    permutations = n,
    strata = metadata$site
  ) |>
    as.data.frame()
} else {
  vegan::adonis2(
    dist ~ habitat * sample_type,
    data = metadata,
    by = "terms",
    permutations = n,
    strata = metadata$site
  ) |>
    as.data.frame()
}

# Inspect output -------------------------------------------------------
print(permanova)

# ======================================================================
# 3. Output
# ======================================================================
write.csv(
  permanova,
  file = here(
    output,
    paste0(data_type, "_", sample_type, "_PERMANOVA_th3.csv")))
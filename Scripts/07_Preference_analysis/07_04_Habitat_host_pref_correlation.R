# 07_07_Habitat_host_pref_correlation.R
#
# Purpose:
#   Test the correlation between habitat preference intensity and host preference 
#   across microbial OTUs using Spearman rank correlation.
#
# Note:
#   Habitat preference is bidirectional (positive/negative values represent different 
#   habitats). The absolute value is used as a measure of preference intensity.
#
# Input:
#   Loaded via "Scripts/07_Preference_analysis/07_00_Setup.R":
#     - data_type: Target data group ("fungi", "prokaryote", etc.)
#     - threshold: Filtering threshold value
#   Files read:
#     - Output/07_Preference_analysis/Proc_data/<data_type>/th<threshold>/preference_data.rds
#
# Output:
#   - Output/07_Preference_analysis/Correlation/<data_type>/th<threshold>/Preference_correlation_result.csv
#
# Analysis:
#   - Spearman rank correlation test between absolute habitat preference and host preference.
#
# R version:
#   R 4.5.3
#
# Packages:
#   here
#   tidyverse
#   stats

# ======================================================================
# 1. Setup & Package Load
# ======================================================================
library(here)
library(tidyverse)
library(stats)

# Load analysis settings -----------------------------------------------
source(here("Scripts", "07_Preference_analysis", "07_00_Setup.R"))

# Define input/output directory paths ----------------------------------
input_file <- here(
  "Output",
  "07_Preference_analysis",
  "Proc_data",
  data_type,
  paste0("th", threshold),
  "preference_data.rds"
)

output <- here(
  "Output",
  "07_Preference_analysis",
  "Correlation",
  data_type,
  paste0("th", threshold)
)

dir.create(output, showWarnings = FALSE, recursive = TRUE)

# ======================================================================
# 2. Load Data
# ======================================================================
preference_df <- readRDS(input_file)

# ======================================================================
# 3. Prepare Data
# ======================================================================
# Habitat preference is bidirectional; absolute values represent intensity.
correlation_data <- preference_df |>
  select(habitat_preference, host_preference) |>
  mutate(habitat_preference = abs(habitat_preference)) |>
  drop_na(habitat_preference, host_preference)

# ======================================================================
# 4. Spearman Rank Correlation
# ======================================================================
# NA removal can leave fewer than two complete pairs.
if (nrow(correlation_data) < 2L) {
  correlation_test <- list(
    estimate = c(rho = NA_real_), statistic = c(S = NA_real_),
    p.value = NA_real_, method = "Spearman's rank correlation rho"
  )
} else {
  correlation_test <- cor.test(
    x = correlation_data$habitat_preference,
    y = correlation_data$host_preference,
    method = "spearman"
  )
}

# Combine test results -------------------------------------------------
correlation_result <- tibble(
  rho     = unname(correlation_test$estimate),
  S       = unname(correlation_test$statistic),
  p.value = correlation_test$p.value,
  method  = correlation_test$method
)

# ======================================================================
# 5. Save Results
# ======================================================================
write.csv(
  correlation_result,
  file = file.path(output, "Preference_correlation_result.csv"),
  row.names = FALSE
)

correlation_result

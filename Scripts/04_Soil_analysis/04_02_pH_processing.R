# 04_02_pH_processing.R
#
# Purpose:
#   Process soil pH data by:
#   (1) calculating descriptive statistics, and
#   (2) testing habitat differences in soil pH.
#
# Input:
#   Data/Soil_analysis/
#     - metadata_patch_level.csv
#
# Output:
#   Output/04_Soil_analysis/02_pH_processing/
#     - pH_summary.csv
#     - pH_test_result.csv
#
# Statistical analysis:
#   - Wilcoxon rank-sum tests were performed separately by site to compare soil pH between habitats.
#   - P-values were adjusted using the Benjamini-Hochberg method.
#
# R version:
#   R 4.5.3
#
# Packages:
#   tidyverse
#   here
#   rstatix

# ======================================================================
# 1. Setup
# ======================================================================
library(tidyverse)
library(here)
library(rstatix)

# Input and output directories ------------------------------------------
input <- here("Data", "Soil_analysis")
output <- here("Output", "04_Soil_analysis", "02_pH_processing")
dir.create(output, showWarnings = FALSE, recursive = TRUE)

# ======================================================================
# 2. Load and prepare pH data
# ======================================================================
# Soil pH measured at the sampling-patch level.
metadata_patch <- read.csv(file.path(input, "metadata_patch_level.csv")) |>
  mutate(habitat = factor(habitat, levels = c("Solfatara field", "Forest edge")))

# ======================================================================
# 3. Descriptive statistics
# ======================================================================
pH_summary <- metadata_patch |>
  group_by(site, habitat) |>
  summarise(
    pH_mean = mean(pH, na.rm = TRUE),
    pH_sd = sd(pH, na.rm = TRUE),
    n = sum(!is.na(pH)),
    .groups = "drop"
  )

write.csv(pH_summary, file.path(output, "pH_summary.csv"), row.names = FALSE)

# ======================================================================
# 4. Wilcoxon rank-sum test
# ======================================================================
# Test habitat differences in soil pH separately by site.
test_pH <- metadata_patch |>
  filter(!is.na(pH)) |>
  group_by(site) |>
  rstatix::wilcox_test(pH ~ habitat) |>
  adjust_pvalue(method = "BH") |>
  add_significance("p.adj")

# ======================================================================
# 5. Calculate Wilcoxon effect sizes
# ======================================================================
effect_df <- metadata_patch |>
  filter(!is.na(pH)) |>
  group_by(site) |>
  wilcox_effsize(pH ~ habitat) |>
  select(site, effsize)

# Combine test results and effect sizes.
test_df <- test_pH |>
  left_join(effect_df, by = "site")

# ======================================================================
# 6. Save statistical results
# ======================================================================
write.csv(test_df, file.path(output, "pH_test_result.csv"), row.names = FALSE)
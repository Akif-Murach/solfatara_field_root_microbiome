# 04_01_Cation_processing.R
#
# Purpose:
#   Process soil exchangeable cation concentrations by:
#   (1) blank correction,
#   (2) conversion to mg kg-1 dry soil,
#   (3) limit-of-detection (LOD) filtering, and
#   (4) statistical comparison between habitats.
#
# Input:
#   Data/Soil_analysis/
#     - KCl+CaCl2_Al_analysis.csv
#     - CH3COONH3_soil_analysis.csv
#     - soil_analysis_metadata.csv
#     - iCAP_PRO_LOD.csv
#     - elements_mode.csv
#
# Output:
#   Output/04_Soil_analysis/01_Cation_processing/
#     - All_ex_cations_corrected_by_blank.csv
#     - LOD_check_result_table.csv
#     - All_ex_cations_long_cor_blank.rds
#     - All_soil_analysis_summary_cor.csv
#     - All_soil_analysis_wilcox_test_result.csv
#
# Notes:
#   - Al was extracted with 1 M KCl + 0.5 M CaCl2.
#   - Other exchangeable cations were extracted with 1 M CH3COONH4 (pH 7).
#   - Concentrations were converted from mg L-1 in the extract to mg kg-1 dry soil.
#   - Two contaminated CH3COONH4 procedural blanks were excluded before calculating the mean blank concentration.
#   - Values below the analytical LOD were converted to NA.
#   - Co was excluded because all samples were below the LOD.
#
# Statistical analysis:
#   - Wilcoxon rank-sum tests were performed separately by site for each element.
#   - P-values were adjusted using the Benjamini-Hochberg method.
#   - Wilcoxon effect sizes were calculated using rstatix.
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
output <- here("Output", "04_Soil_analysis", "01_Cation_processing")
dir.create(output, showWarnings = FALSE, recursive = TRUE)

# ======================================================================
# 2. Load input data
# ======================================================================
# Al analysis: KCl + CaCl2 extraction
Al_data <- read.csv(here(input, "KCl+CaCl2_Al_analysis.csv"))

# Other exchangeable cations: CH3COONH4 extraction
cation_data <- read.csv(here(input, "CH3COONH3_soil_analysis.csv"))

# Soil metadata
metadata <- read.csv(here(input, "soil_analysis_metadata.csv"))

# ======================================================================
# 3. Process Al
# ======================================================================
# 3.1 Blank correction --------------------------------------------------
blank_mean_Al <- Al_data |>
  filter(str_detect(Sample_ID, "blank")) |>
  summarise(mean_Al = mean(Al, na.rm = TRUE)) |>
  pull(mean_Al)

Al_data_cor <- Al_data |>
  filter(!str_detect(Sample_ID, "blank")) |>
  mutate(Al = Al - blank_mean_Al)

# 3.2 Convert to mg kg-1 dry soil --------------------------------------
# Extract concentration: mg L-1
# Extraction: 3 g air-dried soil + 30 mL extractant
# Conversion: mg L-1 × 0.03 L / 0.003 kg = mg kg-1
Al_cal <- Al_data_cor |>
  mutate(Al = Al * 0.03 / 0.003)

# ======================================================================
# 4. Process other exchangeable cations
# ======================================================================
# Select exchangeable cations measured by CH3COONH4 extraction.
# Note: Sample IDs ending in "-a" are standardized before blank filtering.
ex_cation <- cation_data |>
  select(Sample_ID, Na, K, Ca, Mg, Zn, Pb, Ni, Co, Mn) |>
  mutate(Sample_ID = gsub("-a$", "", Sample_ID))

# 4.1 Remove contaminated blanks ---------------------------------------
# Four procedural blanks were prepared.
# Two contaminated blanks (blank1 and blank2) were excluded.
# The remaining blanks were used for blank correction.
ex_cation_fil <- ex_cation |>
  filter(!str_detect(Sample_ID, "blank1|blank2"))

# 4.2 Calculate blank means --------------------------------------------
blank_means <- ex_cation_fil |>
  filter(str_detect(Sample_ID, "blank")) |>
  summarise(across(where(is.numeric), \(x) mean(x, na.rm = TRUE)))

print(blank_means)

# 4.3 Blank correction -------------------------------------------------
ex_cation_cor <- ex_cation_fil |>
  filter(!str_detect(Sample_ID, "blank")) |>
  mutate(across(where(is.numeric), ~ .x - blank_means[[cur_column()]]))

# 4.4 Replace negative concentrations with zero -----------------------
# Negative values after blank correction are set to zero.
ex_cation_cor0 <- ex_cation_cor |>
  mutate(across(where(is.numeric), ~ pmax(.x, 0)))

# 4.5 Convert to mg kg-1 dry soil --------------------------------------
# Extract concentration: mg L-1
# Extraction: 1 g air-dried soil + 20 mL extractant
# Conversion: mg L-1 × 0.02 L / 0.001 kg = mg kg-1
ex_cations_cal <- ex_cation_cor0 |>
  mutate(across(where(is.numeric), ~ .x * 0.02 / 0.001))

# ======================================================================
# 5. Combine Al and other exchangeable cations
# ======================================================================
exchangeable_cations <- ex_cations_cal |>
  left_join(Al_cal, by = "Sample_ID")

write.csv(exchangeable_cations, here(output, "All_ex_cations_corrected_by_blank.csv"), row.names = FALSE)

# ======================================================================
# 6. LOD filtering
# ======================================================================
# Load instrument LOD information and measurement modes.
LOD <- read.csv(here(input, "iCAP_PRO_LOD.csv"))

mode <- read.csv(here(input, "elements_mode.csv")) |>
  mutate(Mode = gsub("Aqueous_", "", Mode))

# 6.1 Reshape LOD table -------------------------------------------------
LOD_long <- LOD |>
  pivot_longer(cols = -Element, names_to = "Mode", values_to = "LOD_Value")

# 6.2 Match elements with their measurement mode ----------------------
# Convert LOD from ug L-1 to mg L-1.
LOD_table <- mode |>
  inner_join(LOD_long, by = c("Element", "Mode")) |>
  mutate(LOD_Value = LOD_Value / 1000) |>
  select(-Mode)

# 6.3 Convert concentration data to long format -----------------------
ex_long <- exchangeable_cations |>
  pivot_longer(cols = -Sample_ID, names_to = "Element", values_to = "Concentration")

# 6.4 Compare concentrations with LOD ---------------------------------
LOD_check <- ex_long |>
  left_join(LOD_table |> select(Element, LOD_Value), by = "Element") |>
  group_by(Element) |>
  mutate(below_LOD = Concentration <= LOD_Value) |>
  ungroup()

# 6.5 Convert values below LOD to NA ----------------------------------
LOD_fil <- LOD_check |>
  mutate(Concentration = if_else(below_LOD, NA_real_, Concentration)) |>
  select(Sample_ID, Element, Concentration) |>
  pivot_wider(names_from = Element, id_cols = Sample_ID, values_from = Concentration)

# 6.6 Summarize LOD detection ------------------------------------------
LOD_summary <- LOD_check |>
  group_by(Element) |>
  summarise(
    Total_Samples = n(),
    Below_LOD_Count = sum(below_LOD, na.rm = TRUE),
    Below_LOD_Rate_Percent = Below_LOD_Count / Total_Samples * 100,
    .groups = "drop"
  )

print(LOD_summary)
write.csv(LOD_summary, here(output, "LOD_check_result_table.csv"), row.names = FALSE)

# ======================================================================
# 7. Exclude elements below LOD in all samples
# ======================================================================
# Co was excluded because 100% of samples were below the LOD.
ex_cation_fil <- LOD_fil |>
  select(-Co)

# ======================================================================
# 8. Prepare data for statistical analysis
# ======================================================================
cation_data_long <- ex_cation_fil |>
  pivot_longer(cols = -Sample_ID, names_to = "Measurement", values_to = "Value") |>
  left_join(metadata, by = "Sample_ID")

saveRDS(cation_data_long, here(output, "All_ex_cations_long_cor_blank.rds"))

# ======================================================================
# 9. Calculate descriptive statistics
# ======================================================================
summary_df <- cation_data_long |>
  group_by(site, Measurement, habitat) |>
  summarise(
    mean = mean(Value, na.rm = TRUE),
    n = sum(!is.na(Value)),
    se = sd(Value, na.rm = TRUE) / sqrt(n),
    .groups = "drop"
  )

write.csv(summary_df, here(output, "All_soil_analysis_summary_cor.csv"), row.names = FALSE)

# ======================================================================
# 10. Wilcoxon rank-sum tests
# ======================================================================
# Habitat differences were tested separately for each site and each element.
test_df <- cation_data_long |>
  filter(!is.na(Value)) |>
  group_by(site, Measurement) |>
  rstatix::wilcox_test(Value ~ habitat) |>
  adjust_pvalue(method = "BH") |>
  mutate(
    p.signif = case_when(
      p.adj <= 0.001 ~ "***",
      p.adj <= 0.01  ~ "**",
      p.adj <= 0.05  ~ "*",
      TRUE           ~ "n.s."))

# ======================================================================
# 11. Calculate Wilcoxon effect sizes
# ======================================================================
effect_df <- cation_data_long |>
  filter(!is.na(Value)) |>
  group_by(site, Measurement) |>
  wilcox_effsize(Value ~ habitat) |>
  select(site, Measurement, effsize)

# Combine test results and effect sizes.
test_df <- test_df |>
  left_join(effect_df, by = c("site", "Measurement"))

# ======================================================================
# 12. Save statistical results
# ======================================================================
write.csv(test_df, here(output, "All_soil_analysis_wilcox_test_result.csv"), row.names = FALSE)
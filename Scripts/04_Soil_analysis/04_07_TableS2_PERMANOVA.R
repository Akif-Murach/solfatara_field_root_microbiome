# 04_07_TableS2_PERMANOVA.R
#
# Purpose:
#   Create a supplementary table of PERMANOVA results for soil chemistry.
#
# Input:
#   Output/04_Soil_analysis/05_Multivariate_analysis/
#     - PERMANOVA_results.csv
#
# Output:
#   Output/04_Soil_analysis/06_Multivariate_analysis_tables/
#     - TableS2_PERMANOVA.docx
#
# R version:
#   R 4.5.3
#
# Packages:
#   dplyr
#   flextable
#   officer
#   here

# ======================================================================
# 1. Setup
# ======================================================================
library(dplyr)
library(flextable)
library(officer)
library(here)

source(here("Function", "make_permanova_table.R"))

# Input and output directories ------------------------------------------
input <- here("Output", "04_Soil_analysis", "05_Multivariate_analysis")
output <- here("Output", "04_Soil_analysis", "06_Multivariate_analysis_tables")
dir.create(output, showWarnings = FALSE, recursive = TRUE)

# ======================================================================
# 2. Load and prepare results
# ======================================================================
results <- read.csv(file.path(input, "PERMANOVA_results.csv"), check.names = FALSE)

table_prep <- results |>
  mutate(Factors = if_else(Factors == "Model", "Habitat", Factors))

# ======================================================================
# 3. Create flextable
# ======================================================================
ft <- make_permanova_table(table_prep)

# ======================================================================
# 4. Create Word document and save
# ======================================================================
doc <- read_docx()

doc <- body_add_par(
  doc,
  "Table S2. Results of permutational multivariate analysis of variance (PERMANOVA) testing for differences in soil chemical composition between habitats.",
  style = "Normal"
)

doc <- body_add_flextable(doc, ft)

print(doc, target = file.path(output, "TableS2_PERMANOVA.docx"))
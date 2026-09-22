# 04_08_TableS3_PERMDISP.R
#
# Purpose:
#   Create a supplementary table of PERMDISP results for soil chemistry.
#
# Input:
#   Output/04_Soil_analysis/05_Multivariate_analysis/
#     - Dispersion_test.csv
#
# Output:
#   Output/04_Soil_analysis/06_Multivariate_analysis_tables/
#     - TableS3_PERMDISP.docx
#
# R version:
#   R 4.5.3
#
# Packages:
#   flextable
#   officer
#   here

# ======================================================================
# 1. Setup
# ======================================================================
library(flextable)
library(officer)
library(here)
library(dplyr)

source(here("Function", "make_permdisp_table.R"))

# Input and output directories ------------------------------------------
input <- here("Output", "04_Soil_analysis", "05_Multivariate_analysis")
output <- here("Output", "04_Soil_analysis", "06_Multivariate_analysis_tables")
dir.create(output, showWarnings = FALSE, recursive = TRUE)

# ======================================================================
# 2. Load results
# ======================================================================
results <- read.csv(file.path(input, "Dispersion_test.csv"), check.names = FALSE)

# ======================================================================
# 3. Create flextable
# ======================================================================
ft <- make_permdisp_table(results, group_name = "Habitat")

# ======================================================================
# 4. Create Word document and save
# ======================================================================
doc <- read_docx()

doc <- body_add_par(
  doc,
  "Table S3. Results of permutational analysis of multivariate dispersions (PERMDISP) testing",
  style = "Normal"
)

doc <- body_add_flextable(doc, ft)

print(doc, target = file.path(output, "TableS3_PERMDISP.docx"))
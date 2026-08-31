# 04_06_TableS1_PCA_loadings.R
#
# Purpose:
#   Create a supplementary table of PCA loadings for soil chemical properties.
#
# Input:
#   Output/04_Soil_analysis/05_Multivariate_analysis/
#     - PCA_loadings.csv
#
# Output:
#   Output/04_Soil_analysis/06_Multivariate_analysis_tables/
#     - TableS1_PCA_Loadings.docx
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

# Input and output directories ------------------------------------------
input <- here("Output", "04_Soil_analysis", "05_Multivariate_analysis")
output <- here("Output", "04_Soil_analysis", "06_Multivariate_analysis_tables")
dir.create(output, showWarnings = FALSE, recursive = TRUE)

# ======================================================================
# 2. Load data
# ======================================================================
loadings <- read.csv(file.path(input, "PCA_loadings.csv"), check.names = FALSE)

# ======================================================================
# 3. Prepare table
# ======================================================================
table_prep <- loadings |>
  mutate(
    Factors = factor(
      Factors,
      levels = c("Na", "K", "Ca", "Mg", "Zn", "Pb", "Ni", "Mn", "Al", "pH")
    )
  ) |>
  mutate(across(where(is.numeric), \(x) round(x, 3)))

# ======================================================================
# 4. Create flextable
# ======================================================================
ft <- flextable(table_prep) |>
  # Format ion names with superscripts
  compose(j = "Factors", i = ~ Factors == "Na", value = as_paragraph("Na", as_sup("+"))) |>
  compose(j = "Factors", i = ~ Factors == "K",  value = as_paragraph("K",  as_sup("+"))) |>
  compose(j = "Factors", i = ~ Factors == "Ca", value = as_paragraph("Ca", as_sup("2+"))) |>
  compose(j = "Factors", i = ~ Factors == "Mg", value = as_paragraph("Mg", as_sup("2+"))) |>
  compose(j = "Factors", i = ~ Factors == "Zn", value = as_paragraph("Zn", as_sup("2+"))) |>
  compose(j = "Factors", i = ~ Factors == "Pb", value = as_paragraph("Pb", as_sup("2+"))) |>
  compose(j = "Factors", i = ~ Factors == "Ni", value = as_paragraph("Ni", as_sup("2+"))) |>
  compose(j = "Factors", i = ~ Factors == "Mn", value = as_paragraph("Mn", as_sup("2+"))) |>
  compose(j = "Factors", i = ~ Factors == "Al", value = as_paragraph("Al", as_sup("3+"))) |>
  # Table formatting
  theme_booktabs() |>
  autofit() |>
  set_table_properties(layout = "autofit", width = 1) |>
  align(align = "center", part = "all") |>
  align(j = 1, align = "left", part = "all") |>
  fontsize(size = 10, part = "all") |>
  bold(part = "header") |>
  font(fontname = "Times New Roman", part = "all")

# ======================================================================
# 5. Create Word document and save
# ======================================================================
caption_bold <- fp_text(font.family = "Times New Roman", font.size = 10, bold = TRUE)
caption_plain <- fp_text(font.family = "Times New Roman", font.size = 10, bold = FALSE)

doc <- read_docx()

doc <- body_add_fpar(
  doc,
  fpar(
    ftext("Table S1. ", prop = caption_bold),
    ftext("Factor loadings of soil chemical properties for the first nine principal components.", prop = caption_plain)
  ),
  style = "Normal"
)

doc <- body_add_flextable(doc, ft)

print(doc, target = file.path(output, "TableS1_PCA_Loadings.docx"))
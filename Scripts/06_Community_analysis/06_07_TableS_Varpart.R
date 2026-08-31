library(dplyr)
library(flextable)
library(officer)
library(here)
#--------------------------------------------------
# 1. Variation partitioning結果の読み込み
#--------------------------------------------------
source(here("Scripts", "06_Community_analysis", "06_00_Setup.R"))

# Input and output directories ------------------------------------------
input <- here("Output", "04_Soil_analysis", "05_Multivariate_analysis")
output <- here("Output", "04_Soil_analysis", "06_Multivariate_analysis_tables")
dir.create(output, showWarnings = FALSE, recursive = TRUE)

#results <- read.csv("th3/Root/varpart_results_summary_th3.csv")|>select(-R2)
results <- read.csv("th3/Root&Soil/varpart_results_summary_th3.csv")|>select(-R2)

#--------------------------------------------------
# 3. flextable
#--------------------------------------------------

ft <- flextable(results) |>
  theme_booktabs() |>
  
  compose(
    part = "header",
    j = "df",
    value = as_paragraph(as_i("df"))
  ) |>
  
  compose(
    part = "header",
    j = "Adj_R2",
    value = as_paragraph("Adjusted ", as_i("R"), as_sup("2"))
  ) |>
  
  autofit() |>
  set_table_properties(layout = "autofit", width = 1) |>
  
  align(align = "center", part = "all") |>
  align(j = "Partition", align = "left", part = "all") |>
  
  fontsize(size = 10, part = "all") |>
  bold(part = "header") |>
  font(fontname = "Times New Roman", part = "all")

#--------------------------------------------------
# 4. Word出力
#--------------------------------------------------

doc <- read_docx()

doc <- body_add_par(
  doc,
  #"Table SX. Results of variation partitioning analysis of root fungal community based on Sorensen dissimilarity.",
  "Table SX. Results of variation partitioning analysis of root and soil fungal community based on Sorensen dissimilarity.",
  style = "Normal"
)

doc <- body_add_flextable(doc, ft)

print(
  doc,
  #target = "th3/Root/ITS_varpart_Table_root.docx"
  target = "th3/Root&Soil/ITS_varpart_Table_root&soil.docx"
)
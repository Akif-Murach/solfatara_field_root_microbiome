# 07_06_Robustness_check.R

# Purpose: Perform robustness check across different filtering thresholds 
# (e.g., th1, th3, th5) using Spearman rank correlation for Z-score matrices (2DP or d').
# All cells with matching row and column names are pooled into one correlation.

# Input: Files read via read_zscore(): Output/07_Preference_analysis/2DP/<data_type>//
# th/2DP_Zvalue_.rds OR Output/07_Preference_analysis/dprime/<data_type>//th/dprime_Zvalue_.rds

# Output: Output/07_Preference_analysis/Robustness_check//<data_type>//
# (* <file_prefix>_th-th_spearman.pdf, * <file_prefix>_cor_testth-th_result.csv)

# R version: R 4.5.3
# Packages: tidyverse, broom, here

# ======================================================================
# 1. Setup & Package Load
# ======================================================================
library(tidyverse)
library(broom)
library(here)

# Settings -------------------------------------------------------------
data_type <- "Prokaryote" # "Prokaryote" or "Fungi"
analysis <- "dprime" # "2DP" or "dprime"
focus <- "host" # "habitat" or "host"
direction <- "host" # "microbe" or "host"
thresholds <- c(1, 3, 5)

output <- here("Output", "07_Preference_analysis", "Robustness_check")
dir.create(output, showWarnings = FALSE, recursive = TRUE)

# ======================================================================
# 2. Functions Definition
# ======================================================================

source(here("Function", "Generate_axis_label.R"))
source(here("Function", "Read_zscore.R"))
source(here("Function", "Prepare_comparison.R"))
source(here("Function", "Run_spearman.R"))
source(here("Function", "Plot_spearman.R"))
source(here("Function", "Run_robustness_check.R"))

# ======================================================================
# 3. Performing Analysis
# ======================================================================
axis_label <- get_axis_label(analysis = analysis, data_type = data_type,
                             focus = focus, direction = direction)
zscore_list <- read_zscore(data_type = data_type, analysis = analysis,
                           focus = focus, direction = direction, 
                           thresholds = thresholds)
threshold_pairs <- combn(thresholds, 2, simplify = FALSE)

condition <- if (analysis == "2DP") focus else direction
analysis_output <- file.path(output, analysis, data_type, condition)
dir.create(analysis_output, showWarnings = FALSE, recursive = TRUE)

file_prefix <- paste0(analysis, "_", data_type, "_", condition)

purrr::walk(threshold_pairs, \(pair) {
  threshold1 <- pair[1]
  threshold2 <- pair[2]
  run_robustness_check(
    mat1 = zscore_list[[paste0("th", threshold1)]],
    mat2 = zscore_list[[paste0("th", threshold2)]],
    threshold1 = threshold1, threshold2 = threshold2,
    label = axis_label, output_dir = analysis_output,
    file_prefix = file_prefix, first_col_only = FALSE,
    zero_lines = (analysis == "2DP")
  )
})

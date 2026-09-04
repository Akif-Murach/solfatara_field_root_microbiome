# 06_02_Variation_partitioning.R
#
# Purpose:
#   Perform variation partitioning to evaluate the unique and shared relative 
#   effects of explanatory factors (Habitat, Host plant, or Sample type) on community structure.
#
# Input:
#   Loaded via "Scripts/06_Community_analysis/06_00_Setup.R":
#     - dist: Distance matrix (or response data)
#     - metadata: Sample metadata
#     - data_type: "Prokaryote" or "Fungi"
#     - sample_type: "Root", "Soil", or "Root&Soil"
#
# Output:
#   Output/03_Community_analysis/02_Variation_partitioning/<data_type>/
#     - <data_type>_<sample_type>_Variation_partitioning_th3.rds
#     - <data_type>_<sample_type>_Variation_partitioning_summary_table_th3.csv
#
# Analysis:
#   - Variation partitioning using vegan::varpart.
#   - Summary table construction with adjusted R-squared, raw R-squared, and degrees of freedom.
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
output <- here("Output", "06_Community_analysis", "02_Variation_partitioning", data_type)
dir.create(output, showWarnings = FALSE, recursive = TRUE)

# ======================================================================
# 2. Variation partitioning calculation
# ======================================================================
vp <- if (sample_type == "Root") {
  vegan::varpart(dist, ~ habitat, ~ host, data = metadata)
} else {
  vegan::varpart(dist, ~ habitat, ~ sample_type, data = metadata)
}

# Save model object ----------------------------------------------------
saveRDS(
  vp,
  file = here(
    output,
    paste0(data_type, "_", sample_type, "_Variation_partitioning_th3.rds")
  )
)

# ======================================================================
# 3. Create and save summary table
# ======================================================================

# Variable labels ------------------------------------------------------
var1 <- "Habitat"
var2 <- if (sample_type == "Root") "Host" else "Sample type"

# Summary table formatting ---------------------------------------------
vp_table <- bind_rows(
  as.data.frame(vp$part$fract),
  as.data.frame(vp$part$indfra)
) |>
  mutate(
    Partition = c(
      rownames(vp$part$fract),
      rownames(vp$part$indfra)
    )
  ) |>
  mutate(
    Partition = recode(
      Partition,
      "[a+c] = X1"      = paste0(var1, " (total)"),
      "[b+c] = X2"      = paste0(var2, " (total)"),
      "[a+b+c] = X1+X2" = paste0(var1, "+", var2),
      "[a] = X1|X2"     = paste0(var1, " (unique)"),
      "[b] = X2|X1"     = paste0(var2, " (unique)"),
      "[c]"             = paste0("Shared (", var1, " \u2229 ", var2, ")"),
      "[d] = Residuals" = "Residual"
    )
  ) |>
  transmute(
    Partition,
    Adj_R2 = sprintf("%.3f", Adj.R.squared),
    R2 = ifelse(is.na(R.squared), "", sprintf("%.3f", R.squared)),
    df = ifelse(is.na(Df), "", Df)
  )
vp_table
# Output CSV -----------------------------------------------------------
write.csv(
  vp_table,
  file = here(
    output,
    paste0(data_type, "_", sample_type, "_Variation_partitioning_summary_table_th3.csv")
  ),
  row.names = FALSE
)

# Read Z-score matrices for all thresholds ----------------------------
read_zscore <- function(data_type, analysis, focus, direction, thresholds) {
  map(thresholds, \(threshold) {
    if (analysis == "2DP") {
      input_dir <- here("Output", "07_Preference_analysis", "2DP",
                        data_type, focus, paste0("th", threshold))
      file <- file.path(input_dir, paste0("2DP_Zvalue_", focus, ".rds"))
    } else if (analysis == "dprime") {
      input_dir <- here("Output", "07_Preference_analysis", "dprime",
                        data_type, direction, paste0("th", threshold))
      file <- file.path(input_dir, paste0("dprime_Zvalue_", direction, ".rds"))
    }
    readRDS(file) |> as.data.frame()
  }) |> set_names(paste0("th", thresholds))
}
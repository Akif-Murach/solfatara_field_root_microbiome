
# Spearman correlation test
run_spearman <- function(df) {
  
  df <- df |>
    filter(
      !is.na(zscore_x),
      !is.na(zscore_y)
    )
  
  cor.test(
    df$zscore_x,
    df$zscore_y,
    method = "spearman"
  )
}


# Spearman correlation test
run_spearman <- function(df) {
  
  df <- df |>
    filter(
      !is.na(zscore_x),
      !is.na(zscore_y)
    )
  
  # Preserve the htest output structure when NA removal leaves too few pairs.
  if (nrow(df) < 2L) {
    return(structure(list(
      statistic = c(S = NA_real_), parameter = NULL,
      p.value = NA_real_, estimate = c(rho = NA_real_),
      null.value = c(rho = 0), alternative = "two.sided",
      method = "Spearman's rank correlation rho",
      data.name = "df$zscore_x and df$zscore_y"
    ), class = "htest"))
  }

  cor.test(
    df$zscore_x,
    df$zscore_y,
    method = "spearman"
  )
}

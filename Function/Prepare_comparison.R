# Prepare two matrices for comparison
prepare_comparison <- function(mat1, mat2, first_col_only = TRUE) {
  common_rows <- intersect(rownames(mat1), rownames(mat2))
  mat1 <- mat1[common_rows, , drop = FALSE]
  mat2 <- mat2[common_rows, , drop = FALSE]
  
  if (first_col_only) {
    return(data.frame(zscore_x = mat1[, 1], zscore_y = mat2[, 1], row.names = common_rows))
  }
  
  mat1 |>
    rownames_to_column("row_label") |>
    pivot_longer(cols = -row_label, names_to = "col_label", values_to = "zscore_x") |>
    unite("combined_name", row_label, col_label, sep = "") |>
    left_join(
      mat2 |>
        rownames_to_column("row_label") |>
        pivot_longer(cols = -row_label, names_to = "col_label", values_to = "zscore_y") |>
        unite("combined_name", row_label, col_label, sep = ""),
      by = "combined_name"
    )
}

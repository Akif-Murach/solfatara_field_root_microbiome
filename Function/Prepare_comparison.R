# Prepare two matrices for comparison
prepare_comparison <- function(mat1, mat2, first_col_only = FALSE) {
  for (mat in list(mat1, mat2)) {
    if (is.null(rownames(mat)) || is.null(colnames(mat)) ||
        anyNA(rownames(mat)) || anyNA(colnames(mat)) ||
        anyDuplicated(rownames(mat)) || anyDuplicated(colnames(mat))) {
      stop("Matrices must have non-missing, unique row and column names.")
    }
  }
  common_rows <- intersect(rownames(mat1), rownames(mat2))
  mat1 <- mat1[common_rows, , drop = FALSE]
  mat2 <- mat2[common_rows, , drop = FALSE]
  
  if (first_col_only) {
    return(data.frame(zscore_x = mat1[, 1], zscore_y = mat2[, 1], row.names = common_rows))
  }
  
  # Align both dimensions by name, then pool all matched cells.
  # Keep identifiers separate to avoid collisions when concatenating names.
  common_cols <- intersect(colnames(mat1), colnames(mat2))
  if (length(common_cols) == 0L) {
    stop("No common column names between thresholds.")
  }
  mat1 <- as.matrix(mat1[, common_cols, drop = FALSE])
  mat2 <- as.matrix(mat2[, common_cols, drop = FALSE])
  data.frame(
    row_label = rep(common_rows, times = length(common_cols)),
    col_label = rep(common_cols, each = length(common_rows)),
    zscore_x = as.vector(mat1),
    zscore_y = as.vector(mat2)
  )
}

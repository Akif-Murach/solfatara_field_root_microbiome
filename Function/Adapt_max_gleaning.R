# Solfatara_Analysis/Function.R
# カスタム関数の定義

#' 重複サンプルからリード数が最大の行を抽出する関数
remove_non_G <- function(data) {
  sample_bases <- gsub("_G_|_G2_|_G3_|_G4_", "_", rownames(data))  
  base_groups <- split(rownames(data), sample_bases)     
  keep_rows <- character()
  
  for (grp in base_groups) {
    if (length(grp) == 1) {
      keep_rows <- c(keep_rows, grp)
    } else {
      sums <- sapply(grp, function(rn) rowSums(data[rn, , drop = FALSE], na.rm = TRUE))
      max_row <- grp[which.max(sums)]  
      keep_rows <- c(keep_rows, max_row)
    }
  }
  return(data[keep_rows, , drop = FALSE])
}


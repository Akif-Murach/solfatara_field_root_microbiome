
# Run one robustness comparison ----------------------------------------
run_robustness_check <- function(mat1, mat2, threshold1, threshold2, label, 
                                 output_dir, file_prefix, first_col_only = TRUE,
                                 zero_lines = TRUE) {
  # Prepare data
  df <- prepare_comparison(mat1 = mat1, mat2 = mat2, first_col_only = first_col_only)
  
  # Spearman correlation
  test <- run_spearman(df)
  print(test)
  
  # Plot
  p <- plot_spearman(
    df = df, test = test,
    x_label = paste0(label, " (threshold = ", threshold1, ")"),
    y_label = paste0(label, " (threshold = ", threshold2, ")"),
    zero_lines = zero_lines
  )
  print(p)
  
  # Output directory creation
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  comparison <- paste0("th", threshold1, "-th", threshold2)
  
  # Save plot
  ggsave(
    filename = file.path(output_dir, paste0(file_prefix, "_", comparison, "_spearman.pdf")),
    plot = p, width = 5, height = 5, units = "in", dpi = 300
  )
  
  # Save correlation test result
  write.csv(tidy(test), file = file.path(output_dir, paste0(file_prefix, "_cor_test", comparison, "_result.csv")), row.names = FALSE)
  invisible(list(data = df, test = test, plot = p))
}


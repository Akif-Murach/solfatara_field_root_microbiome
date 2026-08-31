# Plot Spearman correlation
plot_spearman <- function(df, test, x_label, y_label, zero_lines = TRUE) {
  rho <- unname(test$estimate)
  
  p <- ggplot(df, aes(x = zscore_x, y = zscore_y)) +
    geom_smooth(method = "lm", se = TRUE) +
    geom_point(size = 1.8, alpha = 0.6, color = "black") +
    annotate("text", x = min(df$zscore_x, na.rm = TRUE), y = max(df$zscore_y, na.rm = TRUE),
             hjust = 0, vjust = 1, size = 5, label = paste0("Spearman's rho = ", round(rho, 3))) +
    labs(x = x_label, y = y_label) +
    scale_x_continuous(expand = expansion(mult = 0.05)) +
    scale_y_continuous(expand = expansion(mult = 0.05)) +
    theme_bw(base_size = 15) +
    theme(panel.grid = element_blank(), axis.text = element_text(color = "black"))
  
  if (zero_lines) {
    p <- p +
      geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3, color = "grey50") +
      geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3, color = "grey50")
  }
  return(p)
}

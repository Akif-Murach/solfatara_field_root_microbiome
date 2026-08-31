#' Create a formatted flextable for PERMANOVA results
#'
#' @param results A data frame containing PERMANOVA results.
#'
#' @return A formatted flextable object.
make_permanova_table <- function(results) {
  
  table_prep <- results |>
    rename(
      df = Df,
      SS = SumOfSqs,
      F_val = F,
      P_val = `Pr(>F)`
    ) |>
    mutate(
      across(
        c(SS, R2, F_val),
        \(x) sprintf("%.3f", x)
      ),
      P_val = case_when(
        is.na(P_val) ~ "",
        P_val < 0.001 ~ "< 0.001",
        TRUE ~ sprintf("%.3f", P_val)
      )
    )
  
  flextable(table_prep) |>
    theme_booktabs() |>
    
    compose(
      part = "header",
      j = "df",
      value = as_paragraph(as_i("df"))
    ) |>
    compose(
      part = "header",
      j = "SS",
      value = as_paragraph("SS")
    ) |>
    compose(
      part = "header",
      j = "R2",
      value = as_paragraph(
        as_i("R"),
        as_sup("2")
      )
    ) |>
    compose(
      part = "header",
      j = "F_val",
      value = as_paragraph(
        as_i("F"),
        " value"
      )
    ) |>
    compose(
      part = "header",
      j = "P_val",
      value = as_paragraph(
        as_i("P"),
        "-value"
      )
    ) |>
    
    autofit() |>
    set_table_properties(
      layout = "autofit",
      width = 1
    ) |>
    align(
      align = "center",
      part = "all"
    ) |>
    align(
      j = "Factors",
      align = "left",
      part = "all"
    ) |>
    fontsize(
      size = 10,
      part = "all"
    ) |>
    bold(
      part = "header"
    ) |>
    font(
      fontname = "Times New Roman",
      part = "all"
    )
}
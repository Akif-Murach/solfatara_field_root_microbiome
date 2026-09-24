# make_permdisp_table.R
#
# Purpose:
#   Create a formatted flextable from PERMDISP results.
#
# Input:
#   A data frame containing PERMDISP results.
#
# Output:
#   A formatted flextable object.

make_permdisp_table <- function(
    results,
    group_name = "Habitat",
    p_column = "P_value"
) {
  
  if (!p_column %in% c("P_value", "FDR")) {
    stop(
      "`p_column` must be either 'P_value' or 'FDR'."
    )
  }
  
  table_prep <- results |>
    mutate(
      Factors = recode(
        Factors,
        "Groups" = group_name,
        "Residuals" = "Residuals"
      )
    ) |>
    dplyr::rename(
      df = Df,
      SS = `Sum Sq`,
      MS = `Mean Sq`,
      F_value = F,
      P_value = `Pr(>F)`
    ) |>
    mutate(
      P_report = .data[[p_column]]
    ) |>
    select(
      Factors,
      df,
      SS,
      MS,
      F_value,
      P_report
    ) |>
    mutate(
      across(
        c(SS, MS, F_value),
        \(x) ifelse(
          is.na(x),
          "",
          sprintf("%.3f", x)
        )
      ),
      P_report = case_when(
        is.na(P_report) ~ "",
        P_report < 0.001 ~ "< 0.001",
        P_report < 0.01  ~ "< 0.01",
        P_report < 0.05  ~ "< 0.05",
        TRUE ~ sprintf("%.3f", P_report)
      )
    )
  
  ft <- flextable(table_prep) |>
    theme_booktabs() |>
    compose(
      part = "header",
      j = "df",
      value = as_paragraph(as_i("df"))
    ) |>
    compose(
      part = "header",
      j = "F_value",
      value = as_paragraph(as_i("F"), " value")
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
  
  if (p_column == "FDR") {
    
    ft <- ft |>
      compose(
        part = "header",
        j = "P_report",
        value = as_paragraph("FDR")
      )
    
  } else {
    
    ft <- ft |>
      compose(
        part = "header",
        j = "P_report",
        value = as_paragraph(
          as_i("P"),
          "-value"
        )
      )
  }
  
  ft
}
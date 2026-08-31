# Generate axis label
get_axis_label <- function(analysis, data_type, focus, direction) {
  if (analysis == "2DP") {
    if (focus == "habitat") return("Habitat 2DP z-score")
    if (focus == "host") return("Host 2DP z-score")
  } else if (analysis == "dprime") {
    if (direction == "microbe") return(paste0(data_type, " d' z-score"))
    if (direction == "host") return("Host d' z-score")
  }
  stop("Invalid combination of analysis, focus, or direction.")
}
# 03_03_OTU_histogram.R
library(here)

data_type <- "Prokaryote"  #"Prokaryote" or "Fungi"

# Input
input <- here(
  "Output", "03_Data_filtering", "Seqdata", data_type, "Root",
  paste0(data_type, "_Root_coverage_rared_th1fil.rds"))

seqdata <- readRDS(input)
prevalence <- colSums(seqdata > 0)

# Output
output <- here("Output", "03_Data_filtering", "Figures")
dir.create(output, showWarnings = FALSE, recursive = TRUE)


y_label <- if (data_type == "Prokaryote") {
  "Number of prokaryotic OTUs"
} else if (data_type == "Fungi") {
  "Number of fungal OTUs"
} else {
  paste0("Number of ", data_type, " OTUs")}

pdf(
  file.path(output, paste0(data_type, "_Root_histogram.pdf")),
  width = 8,
  height = 6)

hist(
  prevalence,
  breaks = seq(0, 200, by = 1),
  xlab = "Prevalence (Number of samples)",
  ylab = y_label,
  main="")

dev.off()
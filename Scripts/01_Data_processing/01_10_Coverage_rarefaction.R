# 01_10_Coverage_rarefaction.R
#===========================================================
# Coverage-based rarefaction
#===========================================================
library(here)
library(tidyr)
library(vegan)
library(stringr)
library(parallel)
library(dplyr)

#===========================================================
# Analysis settings
#===========================================================
data_type <- "Prokaryote"   # "Prokaryote" or "Fungi"
sample_type <- "Soil"       # "Root" or "Soil"

# Read-depth threshold for initial filtering
read_threshold <- if (sample_type == "Root") 2000 else 5000

# Input file
input_file <- here("Data", data_type, "Seqdata", "seqOTUtab_filtered.rds")

# Output directory
output_dir <- here("Output", "01_Data_processing", data_type, "Covrfy")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
#===========================================================
# Read data
#===========================================================
OTU_table <- readRDS(input_file)

#===========================================================
# Remain maximum reads and remove others for re-sequenced samples
#===========================================================
source(here("Function", "Adapt_max_gleaning.R"))

rm_seq <- remove_non_G(OTU_table)

# Replace G, G2, G3, and G4 suffixes
rownames(rm_seq) <- gsub("_G4|_G3|_G2|_G", "", rownames(rm_seq))

# Remove negative controls
rmf_seq <- rm_seq[!grepl("nega", rownames(rm_seq)), ]

# Separate root and soil samples
root_seq <- rmf_seq[!grepl("Soil", rownames(rmf_seq)), ]
soil_seq <- rmf_seq[grepl("Soil", rownames(rmf_seq)), ]

# Select sample_type
seq <- if (sample_type == "Root") root_seq else soil_seq

#===========================================================
# Rarefaction curve
#===========================================================

rarecurve_file <- file.path(output_dir,
                            paste0(data_type, "_", sample_type, "_rarecurve.pdf"))

pdf(rarecurve_file, width = 8, height = 6)

rarecurve(seq, label = FALSE)

dev.off()

#===========================================================
# Histogram of sequencing reads for each sample
#===========================================================
sumbdt2 <- rowSums(seq)

hist(sumbdt2, breaks = 100,
  main = paste(data_type, sample_type),
  xlab = "Sequencing reads")

axis(1, at = seq(0, 10000, by = 1000))

#===========================================================
# Filter samples by sequencing depth
#===========================================================

OTU_filtered <- seq[sumbdt2 >= read_threshold, , drop = FALSE]

cat("Number of samples before filtering:", nrow(seq),"\n")
cat("Number of samples after filtering:", nrow(OTU_filtered), "\n")
cat("Read-depth threshold:", read_threshold, "\n")

#===========================================================
# Set up parallel backend
#===========================================================
n_cores <- 9
cl <- makeCluster(n_cores)
clusterExport(cl,
  varlist = c("OTU_filtered"),
  envir = environment())

clusterEvalQ(cl, library(vegan))

#===========================================================
# Compute rarefaction slopes in parallel
#===========================================================
# Calculate rarefaction slopes from 1 to maximum - 1 reads for each sample.
rareslopelist <- parLapply(
  cl,
  1:nrow(OTU_filtered),
  function(i) {
    vegan::rareslope(
      OTU_filtered[i, ],
      1:(sum(OTU_filtered[i, ]) - 1))})

#===========================================================
# Get minimum coverage at maximum sequencing depth
#===========================================================
getmincov <- parSapply(
  cl,
  rareslopelist,
  function(x) x[length(x)])

#===========================================================
# Check coverage threshold
#===========================================================
coverage_threshold <- max(getmincov)
cat("Coverage (%):", (1 - coverage_threshold) * 100, "\n")

#===========================================================
# Get read depths required to reach the coverage threshold
#===========================================================

# Return the minimum sequencing depth satisfying the
# coverage threshold for each rarefaction slope vector.

cvrfun <- function(x) {
  min(which(x <= coverage_threshold)) + 1}

clusterExport(
  cl,
  varlist = c("coverage_threshold", "cvrfun"),
  envir = environment())

cvrrare <- parSapply(
  cl,
  rareslopelist,
  cvrfun)

#===========================================================
# Stop parallel cluster
#===========================================================
stopCluster(cl)

#===========================================================
# Inspect coverage-based sequencing depths
#===========================================================
hist(
  cvrrare,
  main = paste(data_type, sample_type),
  xlab = "Required sequencing depth")

print(cvrrare)

#===========================================================
# Save coverage-based sequencing depths
#===========================================================

coverage_depth_file <- file.path(
  output_dir,
  paste0(data_type, "_", sample_type, "_cov_depth_read.rds"))

saveRDS(cvrrare, coverage_depth_file)

#===========================================================
# Coverage-based rarefaction
#===========================================================
set.seed(1234)
# Randomly subsample each sample according to the sequencing depth
# required to achieve the target coverage.
OTU_covrared <- rrarefy(OTU_filtered, cvrrare)

#===========================================================
# Output
#===========================================================
rarefied_file <- file.path(
  output_dir,paste0(data_type,"_",sample_type,"_coverage_rared.rds"))

saveRDS(OTU_covrared,rarefied_file)

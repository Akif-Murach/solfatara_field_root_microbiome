# 02_01_Process_seq_and_taxa.R
#
# Purpose:
#   Process plant ITS sequences and prepare root and leaf sequences
#   for local BLAST-based plant identification.
#
# Input:
#   Data/Plant/seqOTUtab_filtered.rds
#   Data/Plant/OTUseq_0.97.fasta
#   Function/Adapt_max_gleaning.R
#
# Output:
#   Output/02_Plant_root_identification/Seqdata/processed_root_seqdata.rds
#   Output/02_Plant_root_identification/Seqdata/processed_leaf_seqdata.rds
#   Output/02_Plant_root_identification/Seqdata/Leaf.fasta
#   Output/02_Plant_root_identification/Seqdata/Root.fasta

library(tidyverse)
library(Biostrings)
library(here)

source(here("Function", "Adapt_max_gleaning.R"))
input<-here("Data","Plant","Seqdata")

output<-here("Output","02_Plant_root_identification","Seqdata")
dir.create(output,showWarnings = FALSE,recursive = TRUE)

# ============================================================
# 1. Load sequence data
# ============================================================
data <- readRDS(here(input, "seqOTUtab_filtered.rds")) |>
  as.data.frame()
# OTUs detected in negative controls were checked in the unfiltered
# "seqOTUtab.rds". All were identified as non-Viridiplantae based on
# Claident assignments and manual BLAST searches and were excluded
# by the upstream Viridiplantae filter.

# Remove negative-control samples.
dataf <- data[!grepl("nega", rownames(data)), , drop = FALSE]

# ============================================================
# 2. Remove OTUs absent from biological samples
# ============================================================
# Viridiplantae filtering was completed in 01_08_OTU_clustering.R.
plant_data <- dataf |>
  select(where(~ sum(.) > 0))

# ============================================================
# 3. Select representative samples
# ============================================================
# Keep the sample with the maximum sequencing reads
# for samples with multiple gleaning replicates.
plant_datag <- remove_non_G(plant_data)

sample_id_map <- tibble(
  sequencing_sample_id = rownames(plant_datag),
  Sample_ID = gsub("_(G4|G3|G2|G)_", "_", rownames(plant_datag)))
if (anyDuplicated(sample_id_map$Sample_ID)) {
  stop("Duplicate sample IDs remained after gleaning-ID normalization.")}
rownames(plant_datag) <- sample_id_map$Sample_ID

# ============================================================
# 4. Split root and leaf samples
# ============================================================
is_leaf <- str_detect(rownames(plant_datag), "Leaf")

data_leaf <- plant_datag[is_leaf, ]
data_root <- plant_datag[!is_leaf, ]

data_leaff <- data_leaf[, colSums(data_leaf) > 0]
data_rootf <- data_root[, colSums(data_root) > 0]

# Save processed root and leaf sequence tables.
saveRDS(
  data_leaff,
  here(output, "processed_leaf_seqdata.rds"))

saveRDS(
  data_rootf,
  here(output, "processed_root_seqdata.rds"))

# ============================================================
# 5. Extract root and leaf FASTA files
# ============================================================
fasta <- readDNAStringSet(here(input, "OTUseq_0.97.fasta"))

fasta_leaf <- fasta[names(fasta) %in% colnames(data_leaff)]
fasta_root <- fasta[names(fasta) %in% colnames(data_rootf)]

writeXStringSet(
  fasta_leaf,
  filepath = here(output, "Leaf.fasta"),
  format = "fasta")

writeXStringSet(
  fasta_root,
  filepath = here(output, "Root.fasta"),
  format = "fasta")

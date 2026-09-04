# 02_01_Process_seq_and_taxa.R
#
# Purpose:
#   Process plant ITS sequences and prepare root and leaf sequences
#   for local BLAST-based plant identification.
#
# Input:
#   Data/Plant/seqOTUtab.rds
#   Data/Plant/taxonomy_list.rds
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
data <- readRDS(here(input, "seqOTUtab_filtered.rds")) |> as.data.frame()

# Identify OTUs detected in negative controls.
nega  <- data[grepl("nega", rownames(data)), ]
negaf <- nega[, colSums(nega) > 0]

# OTUs detected in negative controls and negative-control samples were removed,
# because all corresponding OTUs were not Viridiplantae.(Please check "SeqOTUtab.rds")

dataf <- data |>
  select(-any_of(colnames(negaf))) |>
  filter(!grepl("nega",rownames(data)))

# ============================================================
# 2. Keep plant OTUs
# ============================================================
taxa <- readRDS(here(input, "taxonomy_list.rds"))
plant_taxa <- taxa |> dplyr::filter(Kingdom == "Viridiplantae")
plant_data <- dataf |>
  select(any_of(plant_taxa$ID)) |>
  select(where(~ sum(.) > 0))

# ============================================================
# 3. Select representative samples
# ============================================================
# Keep the sample with the maximum sequencing reads
# for samples with multiple gleaning replicates.
plant_datag <- remove_non_G(plant_data)

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

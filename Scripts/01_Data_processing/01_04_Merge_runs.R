# script05_Merge_runs.R
#
# Purpose:
#   Merge denoised sequence tables across sequencing runs
#   and remove chimeric sequences.
#
# Input:
#   Data/<dataset>/Seqdata/<run>/stall_no_rmchimera.rds
#
# Output:
#   Output/01_Data_processing/<dataset>/merge_seqtab.nochim.rds
#
# Datasets:
#   - Prokaryote
#   - Fungi
#   - Plant

# ============================================================
# 0. Packages
# ============================================================
library(dada2)
library(here)

# ============================================================
# 1. General settings
# ============================================================
multithread  <- TRUE
merge_method <- "sum"

# ============================================================
# 2. Dataset-specific configuration
# ============================================================
dataset_config <- list(
  Prokaryote = list(
    sequence_files = c(
      here("Data", "Prokaryote", "Seqdata", "251119", "stall_no_rmchimera.rds"),
      here("Data", "Prokaryote", "Seqdata", "251128", "stall_no_rmchimera.rds"),
      here("Data", "Prokaryote", "Seqdata", "251208", "stall_no_rmchimera.rds"),
      here("Data", "Prokaryote", "Seqdata", "260302", "stall_no_rmchimera.rds")
    )
  ),
  Fungi = list(
    sequence_files = c(
      here("Data", "Fungi", "Seqdata", "241019", "stall_no_rmchimera.rds"),
      here("Data", "Fungi", "Seqdata", "250210", "stall_no_rmchimera.rds"),
      here("Data", "Fungi", "Seqdata", "250218", "stall_no_rmchimera.rds"),
      here("Data", "Fungi", "Seqdata", "250303", "stall_no_rmchimera.rds"),
      here("Data", "Fungi", "Seqdata", "250604", "stall_no_rmchimera.rds"),
      here("Data", "Fungi", "Seqdata", "250703", "stall_no_rmchimera.rds"),
      here("Data", "Fungi", "Seqdata", "250728", "stall_no_rmchimera.rds"),
      here("Data", "Fungi", "Seqdata", "250805", "stall_no_rmchimera.rds"),
      here("Data", "Fungi", "Seqdata", "250828", "stall_no_rmchimera.rds")
    )
  ),
  Plant = list(
    sequence_files = c(
      here("Data", "Plant", "Seqdata", "241019", "stall_no_rmchimera.rds"),
      here("Data", "Plant", "Seqdata", "250210", "stall_no_rmchimera.rds"),
      here("Data", "Plant", "Seqdata", "250218", "stall_no_rmchimera.rds"),
      here("Data", "Plant", "Seqdata", "250728", "stall_no_rmchimera.rds"),
      here("Data", "Plant", "Seqdata", "250805", "stall_no_rmchimera.rds")
    )
  )
)

# ============================================================
# 3. Functions
# ============================================================
load_sequence_tables <- function(paths) {
  tables <- lapply(paths, readRDS)
  names(tables) <- basename(dirname(paths))
  tables
}

merge_sequence_tables <- function(tables) {
  message("Merging sequence tables...")
  merged <- dada2::mergeSequenceTables(tables = tables, repeats = merge_method)
  
  message("Removing chimeras...")
  seqtab_nochim <- dada2::removeBimeraDenovo(
    merged, method = "consensus", multithread = multithread, verbose = TRUE
  )
  seqtab_nochim
}

# ============================================================
# 4. Run
# ============================================================
for (dataset in names(dataset_config)) {
  message("\n========================================")
  message("Processing: ", dataset)
  message("========================================")
  
  output_dir <- here("Output", "01_Data_processing", dataset)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  tables <- load_sequence_tables(dataset_config[[dataset]]$sequence_files)
  seqtab_nochim <- merge_sequence_tables(tables)
  
  message("Total reads after chimera removal: ", sum(seqtab_nochim))
  
  saveRDS(seqtab_nochim, file.path(output_dir, "merge_seqtab.nochim.rds"))
}

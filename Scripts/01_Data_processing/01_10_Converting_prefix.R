# Convert OTU prefixes in the sequence table, taxonomy, and centroid FASTA.
# Prokaryote: X_ -> P_
# Fungi:      X_ -> F_

library(here)
library(Biostrings)

# ===========================================================
# Parameter settings
# ===========================================================

data_type <- "Fungi" # "Prokaryote" or "Fungi"

# Use the table actually deposited and used in downstream analyses.
# Change to "seqOTUtab.rds" only when converting the unfiltered table.
seqtab_file <- "seqOTUtab_filtered.rds"

# Optional command-line argument: directory containing the three target files.
# When omitted, use the original downstream-analysis location under Data.
args <- commandArgs(trailingOnly = TRUE)
input_dir <- if (length(args) >= 1) {
  args[[1]]
} else {
  here("Data", data_type, "Seqdata")
}

seqtab_path <- file.path(input_dir, seqtab_file)
taxonomy_path <- file.path(input_dir, "OTU_merge_taxonomylist.rds")
fasta_path <- file.path(input_dir, "OTUseq_0.97.fasta")

prefix <- switch(
  data_type,
  Prokaryote = "P_",
  Fungi = "F_",
  stop("data_type must be 'Prokaryote' or 'Fungi'.")
)

required_files <- c(seqtab_path, taxonomy_path, fasta_path)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop(
    "Required file(s) not found:\n",
    paste(missing_files, collapse = "\n")
  )
}

# ===========================================================
# Read files and validate OTU IDs
# ===========================================================

seqtab <- readRDS(seqtab_path)
taxonomy <- readRDS(taxonomy_path)
fasta <- Biostrings::readDNAStringSet(fasta_path)

seqtab_ids <- colnames(seqtab)
taxonomy_ids <- rownames(taxonomy)
fasta_ids <- names(fasta)

if (is.null(seqtab_ids)) {
  stop("The sequence table has no column names: ", seqtab_path)
}
if (is.null(taxonomy_ids)) {
  stop("The taxonomy table has no row names: ", taxonomy_path)
}
if (is.null(fasta_ids)) {
  stop("The centroid FASTA has no sequence names: ", fasta_path)
}

if (seqtab_file == "seqOTUtab_filtered.rds" &&
    !setequal(seqtab_ids, taxonomy_ids)) {
  stop("OTU IDs differ between the filtered sequence table and taxonomy table.")
}
if (!all(taxonomy_ids %in% seqtab_ids)) {
  stop("Some taxonomy OTU IDs are absent from the sequence table.")
}
if (!all(seqtab_ids %in% fasta_ids)) {
  stop("Some sequence-table OTU IDs are absent from the centroid FASTA.")
}

# ===========================================================
# Convert prefixes and overwrite the three files
# ===========================================================

all_ids <- c(seqtab_ids, taxonomy_ids, fasta_ids)
target_pattern <- paste0("^", prefix)

if (all(grepl(target_pattern, all_ids))) {
  message("All OTU prefixes have already been converted; no files changed.")
} else {
  if (!all(grepl("^X_", all_ids))) {
    stop("Unexpected or mixed OTU prefixes were found among the three files.")
  }
  
  colnames(seqtab) <- sub("^X_", prefix, seqtab_ids)
  rownames(taxonomy) <- sub("^X_", prefix, taxonomy_ids)
  names(fasta) <- sub("^X_", prefix, fasta_ids)
  
  saveRDS(seqtab, seqtab_path)
  saveRDS(taxonomy, taxonomy_path)
  Biostrings::writeXStringSet(fasta, fasta_path)
  
  message("Converted OTU prefixes for: ", data_type)
}

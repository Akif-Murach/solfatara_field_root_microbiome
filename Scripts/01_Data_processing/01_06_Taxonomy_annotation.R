# 01_06_Taxonomy_annotation.R
#
# Purpose:
#   Prepare consistently named ASV tables and FASTA files for all datasets,
#   then assign Prokaryote and Fungi taxonomy with DADA2. Plant ASVs are
#   assigned later with Claident in 01_07.
#
# Usage:
#   Rscript Scripts/01_Data_processing/01_06_Taxonomy_annotation.R \
#     /path/to/referenceDB
#
# Alternatively, set REFERENCE_DB_DIR before running the script.
#
# Input:
#   Output/01_Data_processing/<dataset>/seqtab_decontam.rds
#
# Output:
#   Output/01_Data_processing/<dataset>/merge_seqtab.rds
#   Output/01_Data_processing/<dataset>/merge_seqtab.txt
#   Output/01_Data_processing/<dataset>/merge_clus_seq.fasta
#   Output/01_Data_processing/{Prokaryote,Fungi}/merge_taxonomylist.rds

library(dada2)
library(seqinr)

# ======================================================================
# 1. Settings
# ======================================================================
args <- commandArgs(trailingOnly = TRUE)
reference_dir <- if (length(args) >= 1) args[[1]] else Sys.getenv(
  "REFERENCE_DB_DIR",
  unset = ""
)

if (!nzchar(reference_dir) || !dir.exists(reference_dir)) {
  stop(
    "Provide the reference database directory as the first argument ",
    "or set REFERENCE_DB_DIR."
  )
}

reference_db <- c(
  Prokaryote = file.path(
    reference_dir,
    "silva_nr99_v138.2_toSpecies_trainset.fa.gz"
  ),
  Fungi = {
    fungi_candidates <- file.path(
      reference_dir,
        "sh_general_release_dynamic_19.02.2025.fasta"
    )
    existing_fungi_reference <- fungi_candidates[file.exists(fungi_candidates)]
    if (length(existing_fungi_reference) > 0) {
      existing_fungi_reference[[1]]
    } else {
      fungi_candidates[[1]]
    }
  }
)

missing_references <- reference_db[!file.exists(reference_db)]
if (length(missing_references) > 0) {
  stop(
    "Reference file(s) not found:\n",
    paste(missing_references, collapse = "\n")
  )
}

multithread <- TRUE
datasets <- c("Prokaryote", "Fungi", "Plant")

# ======================================================================
# 2. Helper functions
# ======================================================================
remove_fungal_prefixes <- function(taxa) {
  prefixes <- c(
    Kingdom = "k__",
    Phylum = "p__",
    Class = "c__",
    Order = "o__",
    Family = "f__",
    Genus = "g__",
    Species = "s__"
  )

  for (rank in intersect(names(prefixes), colnames(taxa))) {
    taxa[, rank] <- sub(paste0("^", prefixes[[rank]]), "", taxa[, rank])
  }

  taxa
}

make_asv_ids <- function(n_asv, dataset) {
  prefix <- switch(dataset, Prokaryote = "P_", Fungi = "F_", Plant = "X_")
  width <- nchar(n_asv)
  sprintf(paste0(prefix, "%0", width, "d"), seq_len(n_asv))
}

# ======================================================================
# 3. Run
# ======================================================================
for (dataset in datasets) {
  message("\nASV preparation and taxonomy: ", dataset)

  output_dir <- file.path("Output", "01_Data_processing", dataset)
  
  seqtab <- readRDS(file.path(output_dir, "seqtab_decontam.rds"))

  taxa <- NULL
  if (dataset %in% names(reference_db)) {
    taxa <- dada2::assignTaxonomy(
      seqtab,
      refFasta = reference_db[[dataset]],
      multithread = multithread
    )
    taxa[is.na(taxa)] <- "Unidentified"

    if (dataset == "Fungi") {
      taxa <- remove_fungal_prefixes(taxa)
    }

    if (!identical(rownames(taxa), colnames(seqtab))) {
      stop("Taxonomy and sequence-table ASV sequences do not match: ", dataset)
    }
  }

  sequences <- colnames(seqtab)
  asv_ids <- make_asv_ids(ncol(seqtab), dataset)
  colnames(seqtab) <- asv_ids
  if (!is.null(taxa)) rownames(taxa) <- asv_ids

  seqinr::write.fasta(
    sequences = as.list(sequences),
    names = asv_ids,
    file.out = file.path(output_dir, "merge_clus_seq.fasta")
  )
  saveRDS(seqtab, file.path(output_dir, "merge_seqtab.rds"))
  write.table(
    cbind(sample = rownames(seqtab), seqtab),
    file.path(output_dir, "merge_seqtab.txt"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  if (!is.null(taxa)) {
    saveRDS(taxa, file.path(output_dir, "merge_taxonomylist.rds"))
  }
}

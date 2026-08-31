# script08_OTU_clustering.R
#
# Purpose:
#   Rename ASVs, cluster ASVs into OTUs at 97% identity,
#   construct OTU tables, and filter OTUs based on taxonomy.
#
# Input:
#   Output/<dataset>/ASV_decontam/seqtab_decontam.rds
#   Output/<dataset>/ASV_decontam/merge_taxonomylist.rds
#
# Output:
#   Output/<dataset>/ASV_decontam/
#     - merge_seqtab.rds
#     - merge_seqtab.txt
#     - merge_clus_seq.fasta
#     - seqOTUtab.rds
#     - seqOTUtab.csv
#     - seqOTUtab_filtered.rds
#     - OTU_merge_taxonomylist.rds
#     - OTUseq_0.97.fasta
#
# OTU clustering:
#   VSEARCH, 97% sequence identity

library(seqinr)
library(here)

otu_identity <- 0.97
vsearch_path <- "vsearch"

dataset_config <- list(
  Prokaryote = list(prefix = "P"),
  Fungi      = list(prefix = "F"),
  Plant      = list(prefix = "Pl")
)

# ============================================================
# 1. Functions
# ============================================================

rename_asvs <- function(seqtab, taxa, prefix) {
  n_asv <- ncol(seqtab)
  width <- nchar(n_asv)
  new_names <- sprintf(paste0(prefix, "_%0", width, "d"), seq_len(n_asv))
  
  fasta_df <- data.frame(sequence = colnames(seqtab), asv = new_names)
  colnames(seqtab) <- new_names
  
  if (!is.null(taxa)) {
    rownames(taxa) <- new_names
  }
  
  list(seqtab = seqtab, taxa = taxa, fasta = fasta_df)
}

write_asv_fasta <- function(fasta_df, output_path) {
  seqinr::write.fasta(
    sequences = as.list(fasta_df[, "sequence"]),
    names     = fasta_df[, "asv"],
    file.out  = output_path
  )
}

cluster_otus <- function(seqtab, output_dir, identity = 0.97) {
  fasta_path     <- file.path(output_dir, "merge_clus_seq.fasta")
  shared_path    <- file.path(output_dir, paste0("ASV_OTU_corestab_", identity, ".txt"))
  centroid_path  <- file.path(output_dir, paste0("OTUseq_", identity, ".fasta"))
  alignment_path <- file.path(output_dir, paste0("seqAlign_", identity, ".txt"))
  
  system2(
    command = vsearch_path,
    args = c(
      "--cluster_fast",     fasta_path,
      "--id",               identity,
      "--mothur_shared_out", shared_path,
      "--centroids",        centroid_path,
      "--msaout",           alignment_path
    )
  )
  
  read.table(shared_path, header = TRUE, row.names = 2)[, -c(1:2), drop = FALSE]
}

make_otu_table <- function(seqtab, otu) {
  otutab <- matrix(
    0, nrow = nrow(seqtab), ncol = ncol(otu),
    dimnames = list(rownames(seqtab), colnames(otu))
  )
  
  for (i in seq_len(ncol(otu))) {
    member_asvs <- rownames(otu)[otu[, i] > 0]
    
    if (length(member_asvs) > 1) {
      otutab[, i] <- rowSums(seqtab[, member_asvs, drop = FALSE])
    } else {
      centroid <- colnames(otu)[i]
      otutab[, i] <- seqtab[, centroid]
    }
  }
  
  otutab
}

filter_otus_by_taxonomy <- function(otutab, taxa, dataset) {
  if (is.null(taxa)) return(otutab)
  
  if (dataset == "Prokaryote") {
    keep <- taxa[, "Family"] != "Mitochondria" & taxa[, "Order"] != "Chloroplast"
  } else if (dataset == "Fungi") {
    keep <- taxa[, "Kingdom"] == "Fungi"
  } else {
    keep <- rep(TRUE, nrow(taxa))
  }
  
  keep[is.na(keep)] <- FALSE
  keep_taxa <- rownames(taxa)[keep]
  
  otutab[, colnames(otutab) %in% keep_taxa, drop = FALSE]
}

# ============================================================
# 2. Run
# ============================================================

for (dataset in names(dataset_config)) {
  message("\n========================================")
  message("OTU clustering: ", dataset)
  message("========================================")
  
  output_dir <- here("Output", dataset, "ASV_decontam")
  seqtab     <- readRDS(file.path(output_dir, "seqtab_decontam.rds"))
  taxa_path  <- file.path(output_dir, "merge_taxonomylist.rds")
  
  taxa <- if (file.exists(taxa_path)) readRDS(taxa_path) else NULL
  
  renamed  <- rename_asvs(seqtab, taxa, dataset_config[[dataset]]$prefix)
  seqtab   <- renamed$seqtab
  taxa     <- renamed$taxa
  fasta_df <- renamed$fasta
  
  write_asv_fasta(fasta_df, file.path(output_dir, "merge_clus_seq.fasta"))
  saveRDS(seqtab, file.path(output_dir, "merge_seqtab.rds"))
  write.table(
    cbind(sample = rownames(seqtab), seqtab),
    file.path(output_dir, "merge_seqtab.txt"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  
  if (!is.null(taxa)) {
    saveRDS(taxa, file.path(output_dir, "merge_taxonomylist.rds"))
  }
  
  otu    <- cluster_otus(seqtab, output_dir, identity = otu_identity)
  otutab <- make_otu_table(seqtab, otu)
  otutab_filtered <- filter_otus_by_taxonomy(otutab, taxa, dataset)
  
  saveRDS(otutab, file.path(output_dir, "seqOTUtab.rds"))
  write.csv(
    cbind(sample = rownames(otutab), otutab),
    file.path(output_dir, "seqOTUtab.csv"),
    row.names = FALSE
  )
  
  if (!is.null(taxa)) {
    otu_taxa <- taxa[rownames(taxa) %in% colnames(otutab_filtered), , drop = FALSE]
    
    saveRDS(otu_taxa, file.path(output_dir, "OTU_merge_taxonomylist.rds"))
    saveRDS(otutab_filtered, file.path(output_dir, "seqOTUtab_filtered.rds"))
  }
}
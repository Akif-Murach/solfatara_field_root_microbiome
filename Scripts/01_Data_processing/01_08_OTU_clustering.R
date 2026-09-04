# 01_08_OTU_clustering.R
#
# Purpose:
#   Cluster ASVs into OTUs at 97% sequence identity, construct OTU tables,
#   extract the taxonomy of representative ASVs, and apply the original
#   dataset-specific taxonomy filters.
#
# Input:
#   Output/01_Data_processing/<dataset>/merge_seqtab.rds
#   Output/01_Data_processing/<dataset>/merge_clus_seq.fasta
#   Output/01_Data_processing/<dataset>/merge_taxonomylist.rds
#
# Output:
#   Output/01_Data_processing/<dataset>/
#     - ASV_OTU_corestab_0.97.txt
#     - OTUseq_0.97.fasta
#     - seqAlign_0.97.txt
#     - seqOTUtab.rds
#     - seqOTUtab.csv
#     - OTU_merge_taxonomylist.rds
#     - seqOTUtab_filtered.rds

library(here)

# ======================================================================
# 1. Settings
# ======================================================================
otu_identity <- 0.97
vsearch_path <- Sys.getenv("VSEARCH_PATH", unset = "vsearch")
datasets <- c("Prokaryote", "Fungi", "Plant")

# ======================================================================
# 2. Helper functions
# ======================================================================
cluster_otus <- function(fasta_path, output_dir, identity) {
  shared_path <- file.path(
    output_dir,
    paste0("ASV_OTU_corestab_", identity, ".txt")
  )
  centroid_path <- file.path(
    output_dir,
    paste0("OTUseq_", identity, ".fasta")
  )
  alignment_path <- file.path(
    output_dir,
    paste0("seqAlign_", identity, ".txt")
  )

  status <- system2(
    command = vsearch_path,
    args = c(
      "--cluster_fast", fasta_path,
      "--id", identity,
      "--mothur_shared_out", shared_path,
      "--centroids", centroid_path,
      "--msaout", alignment_path
    )
  )

  if (!identical(status, 0L)) {
    stop("VSEARCH clustering failed with exit status ", status)
  }

  read.table(
    shared_path,
    header = TRUE,
    row.names = 2,
    check.names = FALSE
  )[, -c(1:2), drop = FALSE]
}

make_otu_table <- function(seqtab, otu_membership) {
  otutab <- matrix(
    0,
    nrow = nrow(seqtab),
    ncol = ncol(otu_membership),
    dimnames = list(rownames(seqtab), colnames(otu_membership))
  )

  for (i in seq_len(ncol(otu_membership))) {
    member_asvs <- rownames(otu_membership)[otu_membership[, i] > 0]

    if (length(member_asvs) > 1) {
      otutab[, i] <- rowSums(seqtab[, member_asvs, drop = FALSE])
    } else {
      centroid <- colnames(otu_membership)[i]
      otutab[, i] <- seqtab[, centroid]
    }
  }

  otutab
}

filter_otus <- function(otutab, taxa, dataset) {
  otu_taxa <- taxa[
    rownames(taxa) %in% colnames(otutab),
    ,
    drop = FALSE
  ]

  keep <- switch(
    dataset,
    Prokaryote = (
      otu_taxa[, "Family"] != "Mitochondria" &
        otu_taxa[, "Order"] != "Chloroplast"
    ),
    Fungi = otu_taxa[, "Kingdom"] == "Fungi",
    Plant = otu_taxa[, "Kingdom"] == "Viridiplantae"
  )
  keep[is.na(keep)] <- FALSE
  otu_taxa <- otu_taxa[keep, , drop = FALSE]

  list(
    taxa = otu_taxa,
    table = otutab[
      ,
      colnames(otutab) %in% rownames(otu_taxa),
      drop = FALSE
    ]
  )
}

# ======================================================================
# 3. Run
# ======================================================================
for (dataset in datasets) {
  message("\nOTU clustering: ", dataset)

  output_dir <- here("Output", "01_Data_processing", dataset)
  seqtab <- readRDS(file.path(output_dir, "merge_seqtab.rds")) |>
    as.matrix()
  taxa <- readRDS(file.path(output_dir, "merge_taxonomylist.rds"))
  fasta_path <- file.path(output_dir, "merge_clus_seq.fasta")

  if (!identical(colnames(seqtab), rownames(taxa))) {
    stop("Taxonomy and ASV table IDs do not match: ", dataset)
  }

  otu_membership <- cluster_otus(
    fasta_path,
    output_dir,
    identity = otu_identity
  )
  otutab <- make_otu_table(seqtab, otu_membership)

  saveRDS(otutab, file.path(output_dir, "seqOTUtab.rds"))
  write.csv(
    cbind(sample = rownames(otutab), otutab),
    file.path(output_dir, "seqOTUtab.csv"),
    row.names = FALSE
  )

  filtered <- filter_otus(otutab, taxa, dataset)
  saveRDS(
    filtered$taxa,
    file.path(output_dir, "OTU_merge_taxonomylist.rds")
  )
  saveRDS(
    filtered$table,
    file.path(output_dir, "seqOTUtab_filtered.rds")
  )
}

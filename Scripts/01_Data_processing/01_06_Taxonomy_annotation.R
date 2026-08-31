# script07_Taxonomy_annotation.R
#
# Purpose:
#   Assign taxonomy to decontaminated ASVs.
#
# Input:
#   Output/<dataset>/ASV_decontam/seqtab_decontam.rds
#
# Output:
#   Output/<dataset>/ASV_decontam/merge_taxonomylist.rds
#
# Notes:
#   - Prokaryote and Fungi taxonomy is assigned using DADA2.
#   - Plant taxonomy is not assigned here because Claident was used.

library(dada2)
library(here)

multithread <- TRUE

reference_db <- list(
  Prokaryote = here("Ref", "silva_nr99_v138.2_toSpecies_trainset.fa.gz"),
  Fungi      = here("Ref", "sh_general_release_dynamic_19.02.2025.fasta")
)

for (dataset in names(reference_db)) {
  message("\n========================================")
  message("Taxonomic annotation: ", dataset)
  message("========================================")
  
  output_dir <- here("Output", dataset, "ASV_decontam")
  seqtab     <- readRDS(file.path(output_dir, "seqtab_decontam.rds"))
  
  taxa <- dada2::assignTaxonomy(
    seqtab,
    refFasta    = reference_db[[dataset]],
    multithread = multithread
  )
  
  taxa[is.na(taxa)] <- "Unidentified"
  
  saveRDS(taxa, file.path(output_dir, "merge_taxonomylist.rds"))
}
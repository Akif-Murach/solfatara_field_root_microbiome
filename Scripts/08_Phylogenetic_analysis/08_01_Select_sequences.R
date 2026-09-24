# 08_01_Select_sequences.R
# Purpose: Select fungal OTUs for phylogenetic analysis based on: (1) availability of preference estimates, and (2) membership in Hyaloscyphaceae.
# Workflow: 1. Load data -> 2. Intersect OTUs -> 3. Extract Hyaloscyphaceae -> 4. Add taxonomy labels -> 5. Combine with ref seqs
# Downstream: MAFFT (MSA) -> trimAl (Trimming) -> MEGA12 (Tree reconstruction)
# Input: Data/Fungi/Seqdata/, Output/04_Preference_analysis/Fungi/Proc_data/
# Output: Output/05_Phylogenetic_analysis/Phylo_data/

# ============================================================
# Packages & Setup
# ============================================================
library(here)
library(Biostrings)
library(tidyverse)

seq_dir <- here("Data", "Fungi", "Seqdata")
pref_dir <- here("Output", "07_Preference_analysis", "Proc_data", "Fungi", "th3")
output_dir <- here("Output", "08_Phylogenetic_analysis", "Phylo_data")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# 1. Load sequence, taxonomy, and preference data
# ============================================================
fasta <- readDNAStringSet(file.path(seq_dir, "OTUseq_0.97.fasta"))
taxa <- readRDS(file.path(seq_dir, "OTU_merge_taxonomylist.rds")) |> 
  as.data.frame()
preference_df <- readRDS(file.path(pref_dir, "preference_data.rds"))

# ============================================================
# 2. Select OTUs available in both datasets
# ============================================================
available_otus <- intersect(preference_df$OTU_ID, names(fasta))
fasta_sub <- fasta[available_otus]
taxa_sub <- taxa[available_otus, , drop = FALSE]

# ============================================================
# 3. Select Hyaloscyphaceae OTUs
# ============================================================
hyaloscyphaceae_otus <- rownames(taxa_sub[taxa_sub$Family == "Hyaloscyphaceae", 
                                          , drop = FALSE])
hyaloscyphaceae_fasta <- fasta_sub[hyaloscyphaceae_otus]

# ============================================================
# 4. Add taxonomy-based labels to FASTA headers
# ============================================================
taxa_hyaloscyphaceae <- taxa_sub[hyaloscyphaceae_otus, , drop = FALSE] |>
  mutate(
    across(c(Genus, Family, Order, Class, Phylum),
           ~ if_else(is.na(.) | grepl("unidentified|Incertae_sedis", ., ignore.case = TRUE), 
                     NA_character_, .)),
    label = coalesce(Genus, Family, Order, Class, Phylum, "Unidentified"))

# Ensure order match & set FASTA headers
taxa_hyaloscyphaceae <- taxa_hyaloscyphaceae[names(hyaloscyphaceae_fasta), , drop = FALSE]
names(hyaloscyphaceae_fasta) <- paste0(taxa_hyaloscyphaceae$label, 
                                       "_[", names(hyaloscyphaceae_fasta), "]")

# ============================================================
# 5. Combine OTU and reference sequences
# ============================================================
reference_fasta <- readDNAStringSet(file.path(seq_dir, 
                                              "Hyaloscyphaceae_refseq.fasta"))
combined_fasta <- c(hyaloscyphaceae_fasta, reference_fasta)

writeXStringSet(combined_fasta, 
                file = file.path(output_dir, "merged_Hyaloscyphaceae_Ref.fasta"))

# ============================================================
# 6. Downstream phylogenetic analysis (External software: MAFFT -> trimAl -> MEGA12)
# ============================================================
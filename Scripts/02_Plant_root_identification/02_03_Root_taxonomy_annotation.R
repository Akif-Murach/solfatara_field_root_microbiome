# 02_03_Root_taxonomy_annotation.R
#
# Purpose:
#   Assign plant identities to root OTUs using local BLAST results
#   against reference leaf OTUs.
#
# Input:
#   Output/02_Plant_root_identification/Seqdata/processed_root_seqdata.rds
#   Output/02_Plant_root_identification/Seqdata/processed_leaf_seqdata.rds
#   Data/Plant/Seqdata/taxonomy_list.rds
#   Data/Plant/Metadata/Leaf_Sample_list.csv
#   Output/02_Plant_root_identification/Local_Blast/Leaf_Root_blast.tsv
#   Output/02_Plant_root_identification/Supplementary_annotation/LBlast_noref_ano.csv
#
# Output:
#   Output/02_Plant_root_identification/local_blast_id97_length100.csv
#   Output/02_Plant_root_identification/Local_blast/Leaf_OTU_identified_list.csv
#   Output/02_Plant_root_identification/Local_blast/Unique_pair_of_Leaf_OTU_and_Genus.csv
#   Output/02_Plant_root_identification/Supplementary_annotation/LBlast_noref.csv
#   Output/02_Plant_root_identification/Seqdata/processed_mergefile.rds

library(tidyverse)
library(here)

input1<-here("Output", "02_Plant_root_identification")
input2<-here("Data", "Plant")
output1<-here(input1,"Local_Blast")
output2<-here(input1,"Supplementary_annotation")
dir.create(output2,showWarnings = FALSE,recursive = TRUE)

# ============================================================
# 1. Load processed sequence data and taxonomy
# ============================================================
data_rootf <- readRDS(here(input1, "Seqdata", "processed_root_seqdata.rds"))
data_leaff <- readRDS(here(input1, "Seqdata", "processed_leaf_seqdata.rds"))
taxa       <- readRDS(here(input2, "Seqdata", "taxonomy_list.rds"))

plant_taxa <- taxa |> dplyr::filter(Kingdom == "Viridiplantae")

# ============================================================
# 2. Process local BLAST results
# ============================================================
blast_colnames <- c(
  "root_id",          # Query sequence ID
  "leaf_id",          # Subject sequence ID
  "pident",           # Percent identity
  "alignment_length", # Alignment length
  "mismatch",         # Number of mismatches
  "gapopen",          # Number of gap openings
  "qstart",           # Query start position
  "qend",             # Query end position
  "sstart",           # Subject start position
  "send",             # Subject end position
  "evalue",           # Expect value
  "bitscore",         # Bit score
  "qlen",             # Query sequence length
  "slen"              # Subject sequence length
)

Local_Blast <- read_tsv(
  here(output1, "Leaf_Root_blast.tsv"),
  col_names = blast_colnames
)

# Retain hits with >=97% sequence identity and an alignment length of >=100 bp.
blast_res_97  <- Local_Blast[Local_Blast$pident >= 97, ]
blast_res_100 <- blast_res_97[blast_res_97$alignment_length >= 100, ]

write.csv(
  blast_res_100,
  here(output1, "local_blast_id97_length100.csv"))

# For each root OTU, select the hit with the highest sequence identity and, when tied, the longest alignment.
result <- blast_res_100 |>
  mutate(across(where(is.list), as.character)) |>
  group_by(root_id) |>
  slice_max(
    order_by = tibble(pident, alignment_length),
    n = 1,
    with_ties = FALSE
  ) |>
  ungroup()

# Extract root-leaf OTU pairs.
LBlast <- result |> select(root_id, leaf_id)

# ============================================================
# 3. Identify leaf OTUs
# ============================================================
Leaf_ident <- read.csv(here(input2, "Metadata", "Leaf_Sample_list.csv"))

# X_001: Rhododendron multiflorum and Rhododendron japonicum -> Rhododendron spp.
# X_004: Eubotryoides grayana (abundant) and Gaultheria adenothrix (rare) -> Eubotryoides grayana

leaf_long <- data_leaff |>
  rownames_to_column(var = "Sample_ID") |>
  pivot_longer(cols = -Sample_ID, names_to = "OTU", values_to = "Count") |>
  filter(Count != 0) |>
  left_join(taxa |> select(ID, Genus), by = c("OTU" = "ID")) |>
  left_join(Leaf_ident |> select(Sample_ID, Identities), by = "Sample_ID")

write.csv(
  leaf_long,
  here(output1, "Leaf_OTU_identified_list.csv"))

# ============================================================
# 4. Assign plant identity to root OTUs
# ============================================================
# X_0004: An identical OTU was detected in Eubotryoides grayana and Gaultheria adenothrix.
# This OTU was treated as Eubotryoides grayana based on the vegetation survey conducted in the field.
refidb <- leaf_long |>
  filter(!str_detect(Identities, "Gaultheria adenothrix")) |>
  mutate(
    Identities_g = str_extract(Identities, "^[^ ]+"),
    Identities_g = str_replace(Identities_g, "Eubotryoides", "Gaultherieae")
  ) |>
  distinct(OTU, Identities_g, .keep_all = TRUE)

write.csv(
  refidb,
  here(output1, "Unique_pair_of_Leaf_OTU_and_Genus.csv"))

refid <- refidb |> select(OTU, Identities_g)

root_identities <- LBlast |>
  left_join(refid, by = c("leaf_id" = "OTU")) |>
  select(root_id, Identities_g) |>
  drop_na(Identities_g) |>
  distinct(root_id, .keep_all = TRUE) |>
  dplyr::rename(OTU = root_id)

# ============================================================
# 5. Extract root OTUs without leaf-reference matches
# ============================================================
noref  <- data_rootf[, !colnames(data_rootf) %in% LBlast$root_id]
noreff <- noref[rowSums(noref) > 0, ]
noref2 <- noreff[, colSums(noreff) > 0]

noref_df <- as.data.frame(ifelse(noref2 > 0, 1, 0))
bdt      <- colSums(noref_df) |> as.data.frame()

taxon <- plant_taxa[plant_taxa$ID %in% rownames(bdt), ]
rownames(taxon) <- taxon$ID

bdt$Genus <- taxon[rownames(bdt), "Genus"]

bdt <- bdt |>
  rownames_to_column(var = "OTU") |>
  dplyr::rename(Count = "colSums(noref_df)")

write.csv(bdt, here(output2, "LBlast_noref.csv"))

# ============================================================
# 6. Annotate remaining OTUs using BLAST top hits
# ============================================================
# If the genus could not be identified from the reference
# leaf samples, the top-hit genus from BLAST was used. 
# (https://blast.ncbi.nlm.nih.gov/)

norefano <- read.csv(here(output2, "LBlast_noref_ano.csv"))

norefano <- norefano |>
  mutate(
    Genus_identities = case_when(
      str_detect(Genus, "unidentified") &
        !str_detect(Top_hit, regex("unidentified", ignore_case = TRUE)) ~ word(Top_hit, 1),
      TRUE ~ Genus
    ),
    Identities_g = case_when(
      str_detect(Genus_identities, "Gaultheria|unidentified_Gaultherieae") ~ "Gaultherieae",
      TRUE ~ Genus_identities))

norefid <- norefano |> select(OTU, Identities_g)

# ============================================================
# 7. Merge plant identity assignments
# ============================================================
mergefile <- rbind(root_identities, norefid)

saveRDS(
  mergefile,
  here(output2, "processed_mergefile.rds"))
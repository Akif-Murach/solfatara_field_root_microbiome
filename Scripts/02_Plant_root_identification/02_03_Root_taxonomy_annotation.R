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
#   Output/02_Plant_root_identification/Local_Blast_id97_length100.csv
#   Output/02_Plant_root_identification/Local_Blast/Leaf_OTU_identified_list.csv
#   Output/02_Plant_root_identification/Local_Blast/Unique_pair_of_Leaf_OTU_and_Genus.csv
#   Output/02_Plant_root_identification/Supplementary_annotation/LBlast_noref.csv
#   Output/02_Plant_root_identification/Supplementary_annotation/processed_mergefile.rds
#   Output/02_Plant_root_identification/Supplementary_annotation/root_leaf_candidates.rds

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
taxa       <- readRDS(here(input2, "Seqdata", "OTU_merge_taxonomylist.rds"))|>
  rownames_to_column(var="OTU")

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
  here(output1, "Local_Blast_id97_length100.csv"))

# Retain ALL qualifying root-leaf pairs. 
# Multiple High-scoring Segment Pairs do not create
# additional candidates or additional read counts.
LBlast <- blast_res_100 |> distinct(root_id, leaf_id)

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
  left_join(taxa |> select(OTU, Genus), by = "OTU") |>
  left_join(Leaf_ident |> select(Sample_ID, Identities), by = "Sample_ID")

write.csv(
  leaf_long,
  here(output1, "Leaf_OTU_identified_list.csv"))

# ============================================================
# 4. Assign plant identity to root OTUs
# ============================================================
# Convert taxonomic assignments to groups eligible for the 90% support calculation.
# Higher-rank assignments are retained in the source columns
# but do not contribute to the numerator, except for Gaultherieae.
normalize_group <- function(x) {
  x <- str_trim(x)
  case_when(
    x %in% c("Eubotryoides", "Gaultheria", "Gaultherieae",
             "unidentified_Gaultherieae") ~ "Gaultherieae",
    is.na(x) | x == "" ~ NA_character_,
    str_detect(x, regex("^(unidentified|unclassified|unknown|uncultured)",
                        ignore_case = TRUE)) ~ NA_character_, TRUE ~ x)}

# Apply the vegetation-based exception only to OTUs shared by both leaf hosts.
shared_gaultherieae <- leaf_long |>
  group_by(OTU) |>
  summarise(shared = all(c("Eubotryoides grayana", "Gaultheria adenothrix") %in%
                          Identities), .groups = "drop") |>
  filter(shared) |> pull(OTU)

refidb <- leaf_long |>
  transmute(
    OTU, reference_host = Identities,
    host_category = case_when(
      Identities %in% c("Rhododendron multiflorum", "Rhododendron japonicum") ~
        "Rhododendron spp.",
      OTU %in% shared_gaultherieae &
        Identities %in% c("Eubotryoides grayana", "Gaultheria adenothrix") ~
        "Eubotryoides grayana",
      TRUE ~ Identities
    ),
    comparison_group = normalize_group(word(Identities, 1))
  ) |> distinct()

write.csv(refidb, here(output1, "Unique_pair_of_Leaf_OTU_and_Genus.csv"),
          row.names = FALSE)

# This audit table may have multiple rows per "root OTU"; NEVER sum reads after
# joining this table directly to the read-abundance table.
leaf_candidates <- LBlast |>
  left_join(refidb, by = c("leaf_id" = "OTU"), relationship = "many-to-many") |>
  dplyr::rename(OTU = root_id) |> distinct()
write.csv(leaf_candidates, here(output2, "Root_leaf_host_candidates.csv"),
          row.names = FALSE)
saveRDS(leaf_candidates, here(output2, "root_leaf_candidates.rds"))

# A process for organizing multiple leaf reference candidates 
# corresponding to each root OTU into a single taxon.
primary <- leaf_candidates |>
  group_by(OTU) |>
  summarise(
    n_groups = n_distinct(comparison_group, na.rm = TRUE),
    comparison_group = if (n_groups == 1L) {
      dplyr::first(comparison_group[!is.na(comparison_group)])
    } else NA_character_,
    .groups = "drop"
  )

# Multiple incompatible leaf groups remain ambiguous. Secondary evidence does
# not silently override that ambiguity. Secondary annotation fills OTUs with
# no usable leaf-library group.
secondary_ids <- setdiff(colnames(data_rootf), primary$OTU[primary$n_groups > 0])
if (anyDuplicated(plant_taxa$ID)) stop("Duplicate OTU IDs in taxonomy_list.rds")
bdt <- tibble(
  OTU = secondary_ids,
  Count = as.integer(colSums(data_rootf[, secondary_ids, drop = FALSE] > 0))
) |> left_join(plant_taxa |> select(OTU, Genus), by = "OTU")
write.csv(bdt, here(output2, "LBlast_noref.csv"), row.names = FALSE)

# Keep the historical, curated annotation file as an explicit input. Join only
# Top_hit: current counts and Claident taxonomy come from the current inputs.
annotation_path <- here(input2, "Supplementary_annotation","LBlast_noref_ano.csv")
manual <- read.csv(annotation_path,stringsAsFactors = FALSE
) |>dplyr::select(OTU, Top_hit)
if (anyDuplicated(manual$OTU)) {
  stop("Duplicate OTUs in LBlast_noref_ano.csv")}
if (!setequal(manual$OTU, bdt$OTU)) {
  stop("OTUs in LBlast_noref_ano.csv do not match ",
       "the current LBlast_noref.csv")}

secondary <- bdt |> left_join(manual, by = "OTU") |>
  mutate(
    claident_group = normalize_group(Genus),
    manual_group = normalize_group(word(Top_hit, 1)))

secondary <- secondary |>
  mutate(
    comparison_group = coalesce(claident_group, manual_group),
    assignment_source = case_when(
      !is.na(claident_group) ~ "claident",
      !is.na(manual_group) ~ "manual_blast",
      TRUE ~ "unresolved"))

# Exactly one row per input OTU: retain unresolved OTUs for the denominator.
mergefile <- tibble(OTU = colnames(data_rootf)) |>
  left_join(primary, by = "OTU") |>
  left_join(secondary |> select(OTU, secondary_group = comparison_group,
                                secondary_source = assignment_source), by = "OTU") |>
  mutate(
    comparison_group = if_else(coalesce(n_groups, 0L) > 1L,
                               NA_character_,
                               coalesce(comparison_group, secondary_group)),
    assignment_source = case_when(
      coalesce(n_groups, 0L) > 1L ~ "ambiguous_leaf",
      coalesce(n_groups, 0L) == 1L ~ "leaf_library",
      TRUE ~ coalesce(secondary_source, "unresolved")
    ),
    comparison_rank = case_when(
      comparison_group == "Gaultherieae" ~ "tribe",
      !is.na(comparison_group) ~ "genus",
      TRUE ~ NA_character_
    )
  ) |> select(OTU, comparison_group, comparison_rank, assignment_source)
stopifnot(!anyDuplicated(mergefile$OTU), nrow(mergefile) == ncol(data_rootf))
write.csv(mergefile, here(output2, "Root_OTU_assignments.csv"), row.names = FALSE)
saveRDS(mergefile, here(output2, "processed_mergefile.rds"))

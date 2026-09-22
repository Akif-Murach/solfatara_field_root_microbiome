# 02_04_Host_filtering.R
#
# Purpose:
#   Filter root samples based on host-assignment confidence and
#   generate host metadata for downstream microbial community analyses.
#
# Input:
#   Output/02_Plant_root_identification/Seqdata/processed_root_seqdata.rds
#   Output/02_Plant_root_identification/Supplementary_annotation/processed_mergefile.rds
#   Output/02_Plant_root_identification/Local_Blast/Unique_pair_of_Leaf_OTU_and_Genus.csv
#   Data/Plant/Metadata/Raw_metadata_sheet.csv
#   Data/Plant/Metadata/Soil_DNA_metadata.csv
#
# Output:
#   Output/02_Plant_root_identification/Metadata/Plant/Host_metadata.csv
#   Output/02_Plant_root_identification/Metadata/Plant/processed_host_metadata.csv
#   Output/02_Plant_root_identification/Metadata/Prokaryote/Metadata/Prokaryote_Root_metadata.csv
#   Output/02_Plant_root_identification/Metadata/Fungi/Fungi_Root_metadata.csv
#   Output/02_Plant_root_identification/Metadata/Prokaryote/Prokaryote_Soil_metadata.csv
#   Output/02_Plant_root_identification/Metadata/Prokaryote/Prokaryote_Root&Soil_metadata.csv
#   Output/02_Plant_root_identification/Metadata/Fungi/Fungi_Soil_metadata.csv
#   Output/02_Plant_root_identification/Metadata/Fungi/Fungi_Root&Soil_metadata.csv

# ============================================================
# Libraries
# ============================================================
library(tidyverse)
library(here)

# ============================================================
# Directories
# ============================================================
input1 <- here("Output", "02_Plant_root_identification")
input2 <- here("Data", "Plant", "Metadata")


output <- here("Output", "02_Plant_root_identification","Metadata")

dir <- list(
  plant = file.path(output, "Plant"),
  prok  = file.path(output, "Prokaryote"),
  fungi = file.path(output, "Fungi"))

lapply(dir,dir.create,showWarnings = FALSE,recursive = TRUE)


# ============================================================
# Functions
# ============================================================
write_meta <- function(data, dir, file) {
  write.csv(data, file.path(dir, file), row.names = FALSE)}

# ============================================================
# Input data
# ============================================================
# Processed root sequencing data
data_rootf <- readRDS(here(input1, "Seqdata", "processed_root_seqdata.rds"))

# Processed OTU-to-host assignment
mergefile  <- readRDS(here(input1, "Supplementary_annotation", "processed_mergefile.rds"))

# Candidate table and one-row-per-OTU assignments are produced together by 02_03.
leaf_candidates <- readRDS(here(input1, "Supplementary_annotation",
                                "root_leaf_candidates.rds"))
stopifnot(!anyDuplicated(mergefile$OTU),
          setequal(mergefile$OTU, colnames(data_rootf)))
if (!"comparison_group" %in% names(mergefile))
  stop("Run the updated 02_03 before 02_04")

#　BEGIN HOST DECISIONS
# Denominator is fixed before any joins and includes unresolved plant reads.
sample_totals <- tibble(Sample_ID = rownames(data_rootf),
                        total_plant_reads = rowSums(data_rootf))
otu_long <- data_rootf |> as.data.frame() |>
  rownames_to_column("Sample_ID") |>
  pivot_longer(-Sample_ID, names_to = "OTU", values_to = "Count") |>
  filter(Count > 0)
annotated_reads <- otu_long |> left_join(mergefile, by = "OTU")
group_support <- annotated_reads |>
  filter(!is.na(comparison_group)) |>
  group_by(Sample_ID, comparison_group) |>
  summarise(group_reads = sum(Count), .groups = "drop") |>
　left_join(sample_totals, by = "Sample_ID") |>
  mutate(prop = group_reads / total_plant_reads)
dominant_groups <- group_support |> filter(prop >= 0.90)
stopifnot(!anyDuplicated(dominant_groups$Sample_ID))

# Candidate rows carry no read counts. Only candidates from OTUs actually
# present in this sample can supply its final leaf-derived host label.
sample_candidates <- otu_long |> select(Sample_ID, OTU) |>
  inner_join(leaf_candidates, by = "OTU", relationship = "many-to-many") |>
  filter(!is.na(comparison_group), !is.na(host_category)) |>
  distinct(Sample_ID, comparison_group, host_category)
host_options <- dominant_groups |>
  inner_join(sample_candidates, by = c("Sample_ID", "comparison_group")) |>
  group_by(Sample_ID) |>
  summarise(
    n_hosts = n_distinct(host_category),
    candidate_hosts = paste(sort(unique(host_category)), collapse = "; "),
    Identities = if (n_hosts == 1L) dplyr::first(host_category) else NA_character_,
    .groups = "drop")

all_decisions <- sample_totals |>
  left_join(dominant_groups |> select(-total_plant_reads), by = "Sample_ID") |>
  left_join(host_options, by = "Sample_ID") |>
  mutate(
    comparison_rank = case_when(
      comparison_group == "Gaultherieae" ~ "tribe",
      !is.na(comparison_group) ~ "genus",
      TRUE ~ NA_character_
    ),
    reason = case_when(
      total_plant_reads <= 0 ~ "no_plant_reads",
      is.na(comparison_group) ~ "below_90_percent",
      is.na(n_hosts) ~ "no_leaf_supported_host",
      n_hosts > 1L ~ "ambiguous_host_category",
      TRUE ~ "retained"))

host_counts <- all_decisions |> filter(reason == "retained") |>
  group_by(Identities) |>
  summarise(n_samples = n_distinct(Sample_ID), .groups = "drop")
all_decisions <- all_decisions |> left_join(host_counts, by = "Identities") |>
  mutate(reason = if_else(reason == "retained" & coalesce(n_samples, 0L) < 10L,
                          "fewer_than_10_samples", reason),
         retained = reason == "retained")
host_info <- all_decisions |> filter(retained)
stopifnot(!anyDuplicated(host_info$Sample_ID))
# END HOST DECISIONS
write_meta(group_support, dir$plant, "Host_group_support.csv")
write_meta(all_decisions, dir$plant, "Host_assignment_decisions.csv")
write_meta(host_info, dir$plant, "Host_metadata.csv")

# ============================================================
# Processed root metadata
# ============================================================
metadata <- read.csv(file.path(input2,"Raw_metadata_sheet.csv"))

host_data <- host_info |>
  select(Sample_ID, Identities) |>
  merge(metadata, by = "Sample_ID") |>
  mutate(
    sample_type = "Root",
    SID = str_extract(Sample_ID, "(?<=_)[^_]+(?=_)") |> 
      str_remove("(?<=[A-Za-z])"))|>
    dplyr::rename(host=Identities)

write_meta(host_data, dir$plant, "processed_host_metadata.csv")

# ============================================================
# Root metadata for downstream microbial analyses
# ============================================================
meta_16s   <- host_data |> mutate(Sample_ID = sub("^Plant", "16S", Sample_ID))
meta_fungi <- host_data |> mutate(Sample_ID = sub("^Plant", "Fungi", Sample_ID))

write_meta(meta_16s, dir$prok, "Prokaryote_Root_metadata.csv")
write_meta(meta_fungi, dir$fungi, "Fungi_Root_metadata.csv")

# ============================================================
# Root + Soil metadata
# ============================================================
soil_info <- read.csv(file.path(input2, "Soil_DNA_metadata.csv")) |>
  mutate(sample_type = "Soil")

root_info <- host_data

# Identify sampling positions represented in both root and soil datasets.
repr <- intersect(soil_info$SID, root_info$SID)

root_info_repr <- root_info |>
  filter(SID %in% repr)

# ============================================================
# Function: Generate Root + Soil metadata
# ============================================================
make_rs_metadata <- function(root, soil, marker, out_dir, prefix) {
  root <- root |>
    mutate(Sample_ID = sub("^Plant", marker, Sample_ID))
  
  soil <- soil |>
    mutate(Sample_ID = paste0("Soil_", marker, "_", SID))
  
  write_meta(soil, out_dir, paste0(prefix, "_Soil_metadata.csv"))
  
  write_meta(
    bind_rows(root, soil),
    out_dir,
    paste0(prefix, "_Root&Soil_metadata.csv")
  )
}

make_rs_metadata(
  root_info_repr,
  soil_info,
  marker = "16S",
  out_dir = dir$prok,
  prefix = "Prokaryote"
)

make_rs_metadata(
  root_info_repr,
  soil_info,
  marker = "Fungi",
  out_dir = dir$fungi,
  prefix = "Fungi"
)

# 02_04_Host_filtering.R
#
# Purpose:
#   Filter root samples based on host-assignment confidence and
#   generate host metadata for downstream microbial community analyses.
#
# Input:
#   Output/02_Plant_root_identification/Seqdata/processed_root_seqdata.rds
#   Output/02_Plant_root_identification/Supplementary_annotation/processed_mergefile.rds
#   Output/02_Plant_root_identification/Local_blast/Unique_pair_of_Leaf_OTU_and_Genus.csv
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

# Local reference database derived from morphologically identified leaf samples
refidb     <- read.csv(here(input1, "Local_blast", "Unique_pair_of_Leaf_OTU_and_Genus.csv"))

# ============================================================
# OTU-level host assignment
# ============================================================
otu_long <- data_rootf |>
  as.data.frame() |>
  rownames_to_column("Sample_ID") |>
  pivot_longer(-Sample_ID, names_to = "OTU", values_to = "Count") |>
  filter(Count > 0) |>
  left_join(mergefile, by = "OTU")

# ============================================================
# Genus-level read abundance
# ============================================================
sample_genus_counts <- otu_long |>
  group_by(Sample_ID, Identities_g) |>
  summarise(genus_reads = sum(Count), .groups = "drop")

sample_genus_prop <- sample_genus_counts |>
  group_by(Sample_ID) |>
  mutate(
    total_reads = sum(genus_reads),
    prop = genus_reads / total_reads
  ) |>ungroup()

# ============================================================
# Host filtering
# ============================================================
# Retain samples in which a single genus accounts for at least 90% of the total reads.
valid_samples <- sample_genus_prop |>
  group_by(Sample_ID) |>
  filter(n() == 1 | max(prop) >= 0.90) |>
  slice_max(order_by = genus_reads, n = 1, with_ties = FALSE) |>
  ungroup() |>
  select(Sample_ID, Identities_g)

filtered <- data_rootf[valid_samples$Sample_ID, ] |>
  rownames_to_column(var = "Sample_ID")

filtered_long <- filtered |>
  pivot_longer(cols = -Sample_ID, names_to = "OTU", values_to = "Count") |>
  filter(Count > 0)

# For each sample, select the OTU with the maximum read count.
filter_max <- filtered_long |>
  group_by(Sample_ID) |>
  slice_max(order_by = Count, n = 1) |>
  ungroup()

# Local reference sequences derived from morphologically identified leaf samples
# were used for host assignment. Claident and manual BLAST assignments were used
# for contamination filtering and taxonomic inspection.
filter_ano <- filter_max |>
  left_join(refidb |> select(OTU, Identities), by = "OTU") |>
  drop_na()

# ============================================================
# Host metadata
# ============================================================
# Retain host identities represented by at least 10 samples.
host_filtering <- filter_ano |>
  group_by(Identities) |>
  count() |>
  filter(n >= 10)

valid_host <- host_filtering$Identities

host_info <- filter_ano |>
  filter(Identities %in% valid_host)

write_meta(host_info, dir$plant, "Host_metadata.csv")

# ============================================================
# Processed root metadata
# ============================================================
metadata <- read.csv(file.path(input2,"Raw_metadata_sheet.csv"))

host_data <- host_info |>
  select(-OTU, -Count) |>
  merge(metadata, by = "Sample_ID") |>
  mutate(
    sample_type = "Root",
    SID = str_extract(Sample_ID, "(?<=_)[^_]+(?=_)") |> 
      str_remove("(?<=[A-Za-z])0"))|>
  rename(host=Identities)

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
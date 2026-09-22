# Normalize sample IDs in raw metadata and three sequence tables.
# ASV/OTU prefixes are assigned in 01_06_Taxonomy_annotation.R.

library(here)

normalize_sample_id <- function(ids) {
  gsub(
    "(^|_)([A-Z]+)0+([1-9][0-9]*)(?=_|$)",
    "\\1\\2\\3",
    ids,perl = TRUE)}

# Plant metadata
metadata_path<- here("Data", "Plant", "Metadata", "Raw_metadata_sheet.csv")
metadata<-read.csv(metadata_path)
metadata$Sample_ID <- normalize_sample_id(metadata$Sample_ID)
write.csv(metadata, metadata_path, row.names = FALSE)

#Soil patch metadata
soil_metapatch_path<-here("Data","Soil_analysis","metadata_patch_level.csv")
soil_metapatch<-read.csv(soil_metapatch_path)
soil_metapatch$Sample_ID<-normalize_sample_id(soil_metapatch$Sample_ID)
write.csv(soil_metapatch, soil_metapatch_path, row.names = FALSE)

# Sequence tables
data_types <- c("Fungi", "Prokaryote", "Plant")

for (data_type in data_types) {
  input_path <- here("Output", "01_Data_processing",
                     data_type, "seqOTUtab_filtered.rds")
  output_dir <- here("Data", data_type, "Seqdata")
  output_path <- file.path(output_dir, "seqOTUtab_filtered.rds")
  
  seqdata <- readRDS(input_path)
  rownames(seqdata) <- normalize_sample_id(rownames(seqdata))
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(seqdata, output_path)}

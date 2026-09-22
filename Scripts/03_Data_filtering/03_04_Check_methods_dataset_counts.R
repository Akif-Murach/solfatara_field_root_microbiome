# check_methods_dataset_counts.R
#
# Purpose:
#   Verify and summarize sample counts, sampling positions (SIDs), and observed 
#   OTU counts for Methods and figure captions at occurrence thresholds 1, 3, and 5.
#   Retain the existing detailed summary restricted to th1 metadata.
#
# Input:
#   Data/<data_type>/Seqdata/
#     - OTU_merge_taxonomylist.rds
#   Output/01_Data_processing/<data_type>/Covrfy/
#     - <data_type>_<sample_type>_coverage_rared.rds
#   Output/02_Plant_root_identification/Metadata/<data_type>/
#     - <data_type>_<sample_type>_metadata.csv
#   Output/03_Data_filtering/Metadata/<data_type>/<sample_type>/
#     - <data_type>_<sample_type>_metadata_th1fil.csv
#
# Output:
#   Output/03_Data_filtering/04_Methods_dataset_counts/
#     - dataset_counts_by_threshold.csv (th1, th3, th5)
#     - dataset_counts.csv (existing th1 summary)
#     - marker_overlap.csv
#     - position_counts.csv
#     - counted_samples.csv
#     - counted_OTUs.csv
#     - excluded_upstream_samples.csv
#
# Analysis & Plot:
#   - Validate sample IDs and sequence/metadata consistency.
#   - Restrict dataset to active samples with >0 reads and filtered th1 metadata.
#   - Calculate taxonomy breakdowns and cross-marker/sample-type position overlaps.
#
# R version:
#   R 4.5.3
#
# Packages:
#   here

# ======================================================================
# 1. Setup Configuration & Output Directory
# ======================================================================
library(here)

output_dir <- here("Output", "03_Data_filtering", "04_Methods_dataset_counts")

data_types <- c("Prokaryote", "Fungi")
sample_types <- c("Root", "Soil")
thresholds <- c(1L, 3L, 5L)

# Check for missing values, empty strings, and duplicate IDs
valid_ids <- function(ids) {
  !is.null(ids) && !anyNA(ids) && all(ids != "") && !anyDuplicated(ids)
}

threshold_summaries <- list()
summaries <- list()
samples <- list()
otus <- list()
exclusions <- list()

# ======================================================================
# 2. Dataset Summarization
# ======================================================================
for (data_type in data_types) {
  taxonomy <- readRDS(here("Data", data_type, "Seqdata", "OTU_merge_taxonomylist.rds"))
  stopifnot(valid_ids(rownames(taxonomy)), "Kingdom" %in% colnames(taxonomy))
  
  for (sample_type in sample_types) {
    dataset <- paste(data_type, sample_type, sep = "_")
    
    # 2-1. Load rarefied OTU matrix and metadata using here()
    seq_file <- here(
      "Output", "01_Data_processing", data_type, "Covrfy",
      paste0(dataset, "_coverage_rared.rds")
    )
    upstream_file <- here(
      "Output", "02_Plant_root_identification", "Metadata",
      data_type, paste0(dataset, "_metadata.csv")
    )
    th1_file <- here(
      "Output", "03_Data_filtering", "Metadata",
      data_type, sample_type, paste0(dataset, "_metadata_th1fil.csv")
    )
    
    seqdata <- readRDS(seq_file)
    upstream_metadata <- read.csv(upstream_file)
    
    # Use th1 metadata for both Root and Soil.
    metadata <- read.csv(th1_file)
    metadata_source <- "03_Data_filtering th1 metadata"
    metadata_samples <- nrow(metadata)
    
    stopifnot(
      is.matrix(seqdata), is.numeric(seqdata),
      all(is.finite(seqdata)), all(seqdata >= 0),
      valid_ids(rownames(seqdata)), valid_ids(colnames(seqdata)),
      valid_ids(metadata$Sample_ID), valid_ids(upstream_metadata$Sample_ID),
      all(c("SID", "site") %in% names(metadata))
    )
    
    # Count each occurrence threshold using the same order as 03_02:
    # match upstream samples -> filter OTUs -> remove empty samples.
    common_ids <- intersect(upstream_metadata$Sample_ID, rownames(seqdata))
    matched <- seqdata[common_ids, , drop = FALSE]
    occurrence <- colSums(matched > 0)

    for (threshold in thresholds) {
      filtered <- matched[, occurrence >= threshold, drop = FALSE]
      filtered <- filtered[rowSums(filtered) > 0, , drop = FALSE]

      # Compare with saved analysis inputs when available. Soil th3/th5
      # can be counted without creating new filtering-stage input files.
      filtered_file <- here(
        "Output", "03_Data_filtering", "Seqdata", data_type, sample_type,
        paste0(dataset, "_coverage_rared_th", threshold, "fil.rds")
      )
      filtered_metadata_file <- here(
        "Output", "03_Data_filtering", "Metadata", data_type, sample_type,
        paste0(dataset, "_metadata_th", threshold, "fil.csv")
      )
      if (file.exists(filtered_file)) {
        saved <- readRDS(filtered_file)
        if (!isTRUE(all.equal(filtered, saved, check.attributes = TRUE))) {
          stop("Recomputed dataset differs from saved input: ", filtered_file)
        }
      }
      if (file.exists(filtered_metadata_file)) {
        saved_metadata <- read.csv(filtered_metadata_file)
        if (!identical(rownames(filtered), saved_metadata$Sample_ID)) {
          stop("Sample IDs/order differ from saved metadata: ", filtered_metadata_file)
        }
      }

      key <- paste0(dataset, "_th", threshold)
      threshold_summaries[[key]] <- data.frame(
        data_type = data_type, sample_type = sample_type,
        threshold = threshold,
        samples = nrow(filtered), OTUs = ncol(filtered)
      )
    }

    # 2-2. Filter for shared samples and exclude rows with zero total reads
    sample_ids <- intersect(metadata$Sample_ID, rownames(seqdata))
    selected <- seqdata[sample_ids, , drop = FALSE]
    selected <- selected[rowSums(selected) > 0, , drop = FALSE]
    sample_ids <- rownames(selected)
    metadata <- metadata[match(sample_ids, metadata$Sample_ID), , drop = FALSE]
    stopifnot(nrow(selected) > 0, !anyNA(metadata$SID), !anyNA(metadata$site))
    
    # 2-3. Count OTUs detected at least once in selected samples by category
    observed_otus <- colnames(selected)[colSums(selected > 0) >= 1]
    if (!all(observed_otus %in% rownames(taxonomy))) {
      stop("Observed OTUs missing taxonomy: ", dataset)
    }
    observed_taxonomy <- taxonomy[observed_otus, , drop = FALSE]
    
    kingdom <- observed_taxonomy[, "Kingdom"]
    unclassified <- is.na(kingdom) |
      kingdom %in% c("", "Unidentified", "Unclassified", "unclassified")
    category <- ifelse(unclassified, "Unclassified", kingdom)
    allowed_categories <- if (data_type == "Prokaryote") {
      c("Bacteria", "Archaea", "Unclassified")
    } else {
      "Fungi"
    }
    if (!all(category %in% allowed_categories)) {
      stop("Unexpected kingdom assignment: ", dataset)
    }
    
    # 2-4. Store summary statistics and corresponding sample ID mapping
    summaries[[dataset]] <- data.frame(
      data_type = data_type, sample_type = sample_type,
      metadata_source = metadata_source, metadata_samples = metadata_samples,
      samples = nrow(selected), positions = length(unique(metadata$SID)),
      OTUs = length(observed_otus),
      Bacteria = sum(category == "Bacteria"), Archaea = sum(category == "Archaea"),
      Unclassified = sum(category == "Unclassified"), Fungi = sum(category == "Fungi")
    )
    
    # Strip marker prefixes to align physical samples while preserving root replicate IDs
    physical_sample <- sub("^Soil_", "", sample_ids)
    physical_sample <- sub("^(16S|Fungi)_", "", physical_sample)
    stopifnot(valid_ids(physical_sample))
    samples[[dataset]] <- data.frame(
      data_type = data_type, sample_type = sample_type,
      Sample_ID = sample_ids, physical_sample = physical_sample,
      SID = metadata$SID, site = metadata$site
    )
    otus[[dataset]] <- data.frame(
      data_type = data_type, sample_type = sample_type,
      OTU = observed_otus, category = category
    )
    
    excluded_ids <- setdiff(upstream_metadata$Sample_ID, sample_ids)
    exclusions[[dataset]] <- data.frame(
      data_type = rep(data_type, length(excluded_ids)),
      sample_type = rep(sample_type, length(excluded_ids)),
      Sample_ID = excluded_ids,
      reason = ifelse(
        !excluded_ids %in% rownames(seqdata),
        "absent from rarefied matrix", "absent from th1 selection or zero reads"
      )
    )
  }
}

# ======================================================================
# 3. Cross-Marker Overlap & Unique Positions Assessment
# ======================================================================
threshold_table <- do.call(rbind, threshold_summaries)
summary_table <- do.call(rbind, summaries)
sample_table <- do.call(rbind, samples)
otu_table <- do.call(rbind, otus)
exclusion_table <- do.call(rbind, exclusions)

overlaps <- list()
for (sample_type in sample_types) {
  prok_samples <- samples[[paste0("Prokaryote_", sample_type)]]$physical_sample
  fungi_samples <- samples[[paste0("Fungi_", sample_type)]]$physical_sample
  overlaps[[sample_type]] <- data.frame(
    sample_type = sample_type,
    both_markers = length(intersect(prok_samples, fungi_samples)),
    either_marker = length(union(prok_samples, fungi_samples)),
    Prokaryote_only = length(setdiff(prok_samples, fungi_samples)),
    Fungi_only = length(setdiff(fungi_samples, prok_samples))
  )
}
overlap_table <- do.call(rbind, overlaps)

# Confirm that each SID uniquely maps to a single site
sites_by_position <- split(sample_table$site, sample_table$SID)
for (sites in sites_by_position) stopifnot(length(unique(sites)) == 1)

root_positions <- unique(sample_table$SID[sample_table$sample_type == "Root"])
soil_positions <- unique(sample_table$SID[sample_table$sample_type == "Soil"])
position_table <- data.frame(
  scope = c("Root", "Soil", "Root_and_Soil"),
  positions = c(
    length(root_positions), length(soil_positions),
    length(union(root_positions, soil_positions))
  )
)

# ======================================================================
# 4. Save Outputs & Display Results
# ======================================================================
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(threshold_table, file.path(output_dir, "dataset_counts_by_threshold.csv"), row.names = FALSE)
write.csv(summary_table, file.path(output_dir, "dataset_counts.csv"), row.names = FALSE)
write.csv(overlap_table, file.path(output_dir, "marker_overlap.csv"), row.names = FALSE)
write.csv(position_table, file.path(output_dir, "position_counts.csv"), row.names = FALSE)
write.csv(sample_table, file.path(output_dir, "counted_samples.csv"), row.names = FALSE)
write.csv(otu_table, file.path(output_dir, "counted_OTUs.csv"), row.names = FALSE)
write.csv(exclusion_table, file.path(output_dir, "excluded_upstream_samples.csv"), row.names = FALSE)

print(threshold_table, row.names = FALSE)
print(summary_table, row.names = FALSE)
print(overlap_table, row.names = FALSE)
print(position_table, row.names = FALSE)
print(exclusion_table, row.names = FALSE)
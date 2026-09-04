# 01_05_Decontamination.R
#
# Purpose:
#   Remove contaminant ASVs from the merged, chimera-filtered Prokaryote
#   and Fungi sequence tables using the prevalence method in decontam.
#   Plant samples are carried forward without decontam, matching the
#   original analysis, because too few negative controls were available.
#
# Input:
#   Output/01_Data_processing/<dataset>/merge_seqtab.nochim.rds
#
# Output:
#   Output/01_Data_processing/<dataset>/seqtab_decontam.rds
#   Output/01_Data_processing/<dataset>/decontam_groups/*.rds
#   Output/01_Data_processing/<dataset>/rm_results/*.csv
#
# The group definitions below reproduce the sample-selection order and
# regular-expression patterns in the original 16S and ITS scripts.

library(decontam)
library(here)

# ======================================================================
# 1. Helper functions
# ======================================================================
run_decontam_group <- function(
    seqtab,
    group_name,
    output_dir,
    samples_to_remove = character()
) {
  is_neg <- grepl("nega", rownames(seqtab))

  if (!any(is_neg)) {
    stop("No negative controls detected in group: ", group_name)
  }
  if (!any(!is_neg)) {
    stop("No biological samples detected in group: ", group_name)
  }

  message(
    "Running decontam: ", group_name,
    " (samples = ", sum(!is_neg),
    ", negative controls = ", sum(is_neg), ")"
  )

  result <- decontam::isContaminant(
    seqtab,
    method = "prevalence",
    neg = is_neg
  )

  removed_asvs <- colnames(seqtab)[result$contaminant]
  neg_samples <- rownames(seqtab)[is_neg]
  real_samples <- rownames(seqtab)[!is_neg]

  removed_summary <- data.frame(
    ASV = removed_asvs,
    prev_neg = colSums(seqtab[neg_samples, removed_asvs, drop = FALSE] > 0),
    prev_sample = colSums(
      seqtab[real_samples, removed_asvs, drop = FALSE] > 0
    ),
    reads_neg = colSums(
      seqtab[neg_samples, removed_asvs, drop = FALSE]
    ),
    reads_sample = colSums(
      seqtab[real_samples, removed_asvs, drop = FALSE]
    )
  )

  write.csv(
    removed_summary,
    file.path(output_dir, "rm_results", paste0(group_name, ".csv")),
    row.names = FALSE
  )

  seqtab_clean <- seqtab[, !result$contaminant, drop = FALSE]
  seqtab_clean <- seqtab_clean[
    !rownames(seqtab_clean) %in% samples_to_remove,
    ,
    drop = FALSE
  ]

  saveRDS(
    seqtab_clean,
    file.path(output_dir, "decontam_groups", paste0(group_name, ".rds"))
  )

  seqtab_clean
}

combine_sequence_tables <- function(seqtabs) {
  all_asvs <- unique(unlist(lapply(seqtabs, colnames)))

  seqtabs_filled <- lapply(seqtabs, function(x) {
    x <- as.data.frame(x)
    missing_asvs <- setdiff(all_asvs, colnames(x))

    if (length(missing_asvs) > 0) {
      x[missing_asvs] <- 0
    }

    x[, all_asvs, drop = FALSE]
  })

  combined <- do.call(rbind, seqtabs_filled)
  as.matrix(combined[, colSums(combined) > 0, drop = FALSE])
}

prepare_output_dirs <- function(dataset) {
  output_dir <- here("Output", "01_Data_processing", dataset)
  dir.create(
    file.path(output_dir, "decontam_groups"),
    recursive = TRUE,
    showWarnings = FALSE
  )
  dir.create(
    file.path(output_dir, "rm_results"),
    recursive = TRUE,
    showWarnings = FALSE
  )
  output_dir
}

# ======================================================================
# 2. Prokaryote decontamination
# ======================================================================
decontam_prokaryote <- function() {
  output_dir <- prepare_output_dirs("Prokaryote")
  seqtab <- readRDS(file.path(output_dir, "merge_seqtab.nochim.rds"))

  root_seq <- seqtab[!grepl("Soil", rownames(seqtab)), , drop = FALSE]

  root1 <- root_seq[!grepl("G", rownames(root_seq)), , drop = FALSE]
  extraction_negatives <- rownames(root1)[
    grepl("nega", rownames(root1)) & !grepl("PCR", rownames(root1))
  ]
  root1_clean <- run_decontam_group(
    root1, "decontam_rmASVs_root1", output_dir
  )

  root2 <- root_seq[
    grepl("G|pheno|ctab|freeze", rownames(root_seq)),
    ,
    drop = FALSE
  ]
  root2_clean <- run_decontam_group(
    root2,
    "decontam_rmASVs_root2",
    output_dir,
    samples_to_remove = extraction_negatives
  )

  soil_seq <- seqtab[grepl("Soil", rownames(seqtab)), , drop = FALSE]
  soil1 <- soil_seq[!grepl("G", rownames(soil_seq)), , drop = FALSE]
  soil1_clean <- run_decontam_group(
    soil1, "decontam_rmASVs_soil1", output_dir
  )

  soil2_pool <- seqtab[
    grepl("G|ex", rownames(seqtab)),
    ,
    drop = FALSE
  ]
  soil2 <- soil2_pool[
    grepl("Soil|nega", rownames(soil2_pool)),
    ,
    drop = FALSE
  ]
  soil2_clean <- run_decontam_group(
    soil2, "decontam_rmASVs_soil2", output_dir
  )

  combine_sequence_tables(
    list(root1_clean, root2_clean, soil1_clean, soil2_clean)
  )
}

# ======================================================================
# 3. Fungi decontamination
# ======================================================================
decontam_fungi <- function() {
  output_dir <- prepare_output_dirs("Fungi")
  seqtab <- readRDS(file.path(output_dir, "merge_seqtab.nochim.rds"))

  root_seq <- seqtab[!grepl("Soil", rownames(seqtab)), , drop = FALSE]

  root1 <- root_seq[!grepl("G", rownames(root_seq)), , drop = FALSE]
  extraction_negatives <- rownames(root1)[
    grepl("nega", rownames(root1)) & !grepl("PCR", rownames(root1))
  ]
  root1_clean <- run_decontam_group(
    root1, "decontam_rmASVs_root1", output_dir
  )

  root_patterns <- c(
    root2 = "G_|pheno|ctab|freeze",
    root3 = "G2|pheno|ctab|freeze",
    root4 = "G3|pheno|ctab|freeze",
    root5 = "G4|pheno|ctab|freeze"
  )

  other_roots <- lapply(names(root_patterns), function(group_name) {
    group_seqtab <- root_seq[
      grepl(root_patterns[[group_name]], rownames(root_seq)),
      ,
      drop = FALSE
    ]

    run_decontam_group(
      group_seqtab,
      paste0("decontam_rmASVs_", group_name),
      output_dir,
      samples_to_remove = extraction_negatives
    )
  })

  soil <- seqtab[grepl("Soil", rownames(seqtab)), , drop = FALSE]
  soil_clean <- run_decontam_group(
    soil, "decontam_rmASVs_soil", output_dir
  )

  combine_sequence_tables(c(list(root1_clean), other_roots, list(soil_clean)))
}

# ======================================================================
# 4. Run
# ======================================================================
for (dataset in c("Prokaryote", "Fungi")) {
  message("\nDecontamination: ", dataset)
  output_dir <- prepare_output_dirs(dataset)

  seqtab_decontam <- switch(
    dataset,
    Prokaryote = decontam_prokaryote(),
    Fungi = decontam_fungi()
  )

  saveRDS(seqtab_decontam, file.path(output_dir, "seqtab_decontam.rds"))
}

# Plant was not subjected to decontam in the original analysis.
plant_dir <- here("Output", "01_Data_processing", "Plant")
plant_seqtab <- readRDS(file.path(plant_dir, "merge_seqtab.nochim.rds"))
saveRDS(plant_seqtab, file.path(plant_dir, "seqtab_decontam.rds"))
message("Plant: decontam skipped; chimera-filtered table carried forward.")

# Remove contaminant ASVs from Prokaryote and Fungi using decontam prevalence.
# Input:  Output/01_Data_processing/<dataset>/merge_seqtab.nochim.rds
# Output: seqtab_decontam.rds, decontam_groups/*.rds, rm_results/*.csv,
#         excluded_negative_samples.csv (and a duplicate report on error).
# Preserve the original group definitions and decontam defaults.
# Negative controls are used for detection, then removed before combination.
# Plant is carried forward unchanged.

library(here)
library(decontam)

negative_pattern <- "nega"

is_negative_sample <- function(sample_ids) {
  grepl(negative_pattern, sample_ids)
}

assert_unique_sample_ids <- function(seqtab, context = "sequence table") {
  sample_ids <- rownames(seqtab)
  if (is.null(sample_ids)) {
    stop("Sample IDs (rownames) are missing in: ", context)
  }

  duplicated_ids <- unique(sample_ids[duplicated(sample_ids)])
  if (length(duplicated_ids) > 0) {
    stop("Duplicated sample IDs detected in ", context, ":\n", paste(duplicated_ids, collapse = "\n"))
  }
  invisible(TRUE)
}

run_decontam_group <- function(seqtab, group_name, output_dir, samples_to_remove = character()) {
  assert_unique_sample_ids(seqtab, context = paste0("decontam group '", group_name, "'"))

  is_neg <- is_negative_sample(rownames(seqtab))
  if (!any(is_neg)) stop("No negative controls detected in group: ", group_name)
  if (!any(!is_neg)) stop("No biological samples detected in group: ", group_name)

  result <- decontam::isContaminant(seqtab, method = "prevalence", neg = is_neg)

  removed_asvs <- colnames(seqtab)[result$contaminant]
  neg_samples  <- rownames(seqtab)[is_neg]
  real_samples <- rownames(seqtab)[!is_neg]

  removed_summary <- data.frame(
    ASV          = removed_asvs,
    prev_neg     = colSums(seqtab[neg_samples,  removed_asvs, drop = FALSE] > 0),
    prev_sample  = colSums(seqtab[real_samples, removed_asvs, drop = FALSE] > 0),
    reads_neg    = colSums(seqtab[neg_samples,  removed_asvs, drop = FALSE]),
    reads_sample = colSums(seqtab[real_samples, removed_asvs, drop = FALSE])
  )

  write.csv(
    removed_summary,
    file.path(output_dir, "rm_results", paste0(group_name, ".csv")),
    row.names = FALSE
  )

  seqtab_clean <- seqtab[, !result$contaminant, drop = FALSE]

  # Remove selected extraction controls only after contaminant detection.
  seqtab_clean <- seqtab_clean[!rownames(seqtab_clean) %in% samples_to_remove, , drop = FALSE]

  saveRDS(seqtab_clean, file.path(output_dir, "decontam_groups", paste0(group_name, ".rds")))
  seqtab_clean
}

combine_sequence_tables <- function(seqtabs, output_dir, dataset) {
  if (is.null(names(seqtabs)) || any(names(seqtabs) == "")) {
    stop("All sequence tables supplied to combine_sequence_tables() must have group names.")
  }

  negative_records <- do.call(
    rbind,
    lapply(names(seqtabs), function(group_name) {
      x <- seqtabs[[group_name]]
      neg_ids <- rownames(x)[is_negative_sample(rownames(x))]
      if (length(neg_ids) == 0) return(NULL)

      data.frame(
        dataset = dataset,
        group   = group_name,
        sample  = neg_ids,
        stringsAsFactors = FALSE
      )
    })
  )

  if (!is.null(negative_records) && nrow(negative_records) > 0) {
    write.csv(negative_records, file.path(output_dir, "excluded_negative_samples.csv"), row.names = FALSE)
  }

  # Shared negative controls must be excluded before combining groups.
  seqtabs_biological <- lapply(seqtabs, function(x) {
    x[!is_negative_sample(rownames(x)), , drop = FALSE]
  })

  sample_membership <- do.call(
    rbind,
    lapply(names(seqtabs_biological), function(group_name) {
      x <- seqtabs_biological[[group_name]]
      if (nrow(x) == 0) return(NULL)

      data.frame(
        sample = rownames(x),
        group  = group_name,
        stringsAsFactors = FALSE
      )
    })
  )

  duplicated_biological_ids <- unique(
    sample_membership$sample[
      duplicated(sample_membership$sample)
    ]
  )

  if (length(duplicated_biological_ids) > 0) {
    duplicate_report <- sample_membership[sample_membership$sample %in% duplicated_biological_ids, , drop = FALSE]
    duplicate_report <- duplicate_report[order(duplicate_report$sample, duplicate_report$group), , drop = FALSE]

    write.csv(duplicate_report, file.path(output_dir, "ERROR_duplicated_biological_samples.csv"), row.names = FALSE)

    duplicate_text <- vapply(
      split(duplicate_report$group, duplicate_report$sample),
      function(groups) paste(unique(groups), collapse = ", "),
      character(1)
    )

    stop(
      dataset, ": duplicated biological sample IDs remain after negative-control removal. ",
      "No samples were renamed and no final table was written.\n",
      paste(paste0(names(duplicate_text), " -> ", duplicate_text), collapse = "\n"),
      "\nSee: ", file.path(output_dir, "ERROR_duplicated_biological_samples.csv")
    )
  }

  # Align ASV columns and fill missing columns with zero.
  all_asvs <- unique(unlist(lapply(seqtabs_biological, colnames), use.names = FALSE))

  seqtabs_filled <- lapply(seqtabs_biological, function(x) {
    x <- as.matrix(x)
    filled <- matrix(0, nrow = nrow(x), ncol = length(all_asvs), dimnames = list(rownames(x), all_asvs))
    if (nrow(x) > 0 && ncol(x) > 0) {
      filled[, colnames(x)] <- x
    }
    filled
  })

  combined <- do.call(rbind, seqtabs_filled)

  combined <- combined[, colSums(combined) > 0, drop = FALSE]

  combined
}

prepare_output_dirs <- function(dataset) {
  output_dir <- here("Output", "01_Data_processing", dataset)
  dir.create(file.path(output_dir, "decontam_groups"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(output_dir, "rm_results"), recursive = TRUE, showWarnings = FALSE)
  output_dir
}

decontam_prokaryote <- function(output_dir) {

  seqtab <- readRDS(file.path(output_dir, "merge_seqtab.nochim.rds"))
  assert_unique_sample_ids(seqtab, context = "Prokaryote merge_seqtab.nochim.rds")

  root_seq <- seqtab[!grepl("Soil", rownames(seqtab)), , drop = FALSE]
  root1    <- root_seq[!grepl("G", rownames(root_seq)), , drop = FALSE]

  extraction_negatives <- rownames(root1)[
    is_negative_sample(rownames(root1)) & !grepl("PCR", rownames(root1))
  ]

  root1_clean <- run_decontam_group(root1, "decontam_rmASVs_root1", output_dir)

  root2 <- root_seq[grepl("G|pheno|ctab|freeze", rownames(root_seq)), , drop = FALSE]
  root2_clean <- run_decontam_group(
    root2, "decontam_rmASVs_root2", output_dir, samples_to_remove = extraction_negatives
  )

  soil_seq    <- seqtab[grepl("Soil", rownames(seqtab)), , drop = FALSE]
  soil1       <- soil_seq[!grepl("G", rownames(soil_seq)), , drop = FALSE]
  soil1_clean <- run_decontam_group(soil1, "decontam_rmASVs_soil1", output_dir)

  soil2_pool  <- seqtab[grepl("G|ex", rownames(seqtab)), , drop = FALSE]
  soil2       <- soil2_pool[grepl("Soil|nega", rownames(soil2_pool)), , drop = FALSE]
  soil2_clean <- run_decontam_group(soil2, "decontam_rmASVs_soil2", output_dir)

  combine_sequence_tables(
    seqtabs = list(
      root1 = root1_clean,
      root2 = root2_clean,
      soil1 = soil1_clean,
      soil2 = soil2_clean
    ),
    output_dir = output_dir,
    dataset    = "Prokaryote"
  )
}

decontam_fungi <- function(output_dir) {

  seqtab <- readRDS(file.path(output_dir, "merge_seqtab.nochim.rds"))
  assert_unique_sample_ids(seqtab, context = "Fungi merge_seqtab.nochim.rds")

  root_seq <- seqtab[!grepl("Soil", rownames(seqtab)), , drop = FALSE]

  root1 <- root_seq[!grepl("G", rownames(root_seq)), , drop = FALSE]

  extraction_negatives <- rownames(root1)[
    is_negative_sample(rownames(root1)) & !grepl("PCR", rownames(root1))
  ]

  root1_clean <- run_decontam_group(root1, "decontam_rmASVs_root1", output_dir)

  root_patterns <- c(
    root2 = "G_|pheno|ctab|freeze",
    root3 = "G2|pheno|ctab|freeze",
    root4 = "G3|pheno|ctab|freeze",
    root5 = "G4|pheno|ctab|freeze"
  )

  other_roots <- lapply(names(root_patterns), function(group_name) {
    group_seqtab <- root_seq[grepl(root_patterns[[group_name]], rownames(root_seq)), , drop = FALSE]
    run_decontam_group(
      group_seqtab,
      paste0("decontam_rmASVs_", group_name),
      output_dir,
      samples_to_remove = extraction_negatives
    )
  })
  names(other_roots) <- names(root_patterns)

  soil       <- seqtab[grepl("Soil", rownames(seqtab)), , drop = FALSE]
  soil_clean <- run_decontam_group(soil, "decontam_rmASVs_soil", output_dir)

  combine_sequence_tables(
    seqtabs    = c(list(root1 = root1_clean), other_roots, list(soil = soil_clean)),
    output_dir = output_dir,
    dataset    = "Fungi"
  )
}

for (dataset in c("Prokaryote", "Fungi")) {

  output_dir <- prepare_output_dirs(dataset)

  seqtab_decontam <- switch(
    dataset,
    Prokaryote = decontam_prokaryote(output_dir),
    Fungi      = decontam_fungi(output_dir)
  )

  saveRDS(seqtab_decontam, file.path(output_dir, "seqtab_decontam.rds"))
}

# Plant: carry forward unchanged, including any negative controls.
plant_dir    <- here("Output", "01_Data_processing", "Plant")
plant_seqtab <- readRDS(file.path(plant_dir, "merge_seqtab.nochim.rds"))

assert_unique_sample_ids(plant_seqtab, context = "Plant merge_seqtab.nochim.rds")
saveRDS(plant_seqtab, file.path(plant_dir, "seqtab_decontam.rds"))


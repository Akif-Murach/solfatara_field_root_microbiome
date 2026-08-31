# script06_Decontamination.R
#
# Purpose:
#   Remove contaminant ASVs using decontam.
#
# Input:
#   Output/<dataset>/ASV_decontam/merge_seqtab.nochim.rds
#
# Output:
#   Output/<dataset>/ASV_decontam/seqtab_decontam.rds
#   Output/<dataset>/ASV_decontam/rm_results/
#
# Notes:
#   - Decontam is applied to Prokaryote and Fungi.
#   - Plant ITS is not subjected to decontam because
#     negative controls were insufficient.

library(decontam)
library(here)

# ============================================================
# 1. Functions
# ============================================================

run_decontam <- function(
    seqtab,
    neg_pattern = "nega",
    remove_pattern = NULL,
    group_name = NULL,
    output_dir = NULL) {
  
  sample_names <- rownames(seqtab)
  is_neg <- grepl(neg_pattern, sample_names)
  
  if (sum(is_neg) == 0) {
    stop("No negative controls detected in group: ", group_name)
  }
  if (sum(!is_neg) == 0) {
    stop("No biological samples detected in group: ", group_name)
  }
  
  metadata <- data.frame(is.neg = is_neg, row.names = sample_names)
  
  message(
    "Running decontam: ", group_name,
    " (samples = ", sum(!is_neg),
    ", negative controls = ", sum(is_neg), ")"
  )
  
  contam_result <- decontam::isContaminant(
    seqtab,
    method = "prevalence",
    neg = metadata$is.neg
  )
  
  seqtab_clean <- seqtab[, !contam_result$contaminant, drop = FALSE]
  
  if (!is.null(remove_pattern)) {
    remove_samples <- grepl(remove_pattern, rownames(seqtab_clean))
    seqtab_clean <- seqtab_clean[!remove_samples, , drop = FALSE]
  }
  
  if (!is.null(output_dir)) {
    removed_asvs <- colnames(seqtab)[contam_result$contaminant]
    
    if (length(removed_asvs) > 0) {
      neg_samples  <- rownames(seqtab)[is_neg]
      real_samples <- rownames(seqtab)[!is_neg]
      
      removed_summary <- data.frame(
        ASV          = removed_asvs,
        prev_neg     = colSums(seqtab[neg_samples, removed_asvs, drop = FALSE] > 0),
        prev_sample  = colSums(seqtab[real_samples, removed_asvs, drop = FALSE] > 0),
        reads_neg    = colSums(seqtab[neg_samples, removed_asvs, drop = FALSE]),
        reads_sample = colSums(seqtab[real_samples, removed_asvs, drop = FALSE])
      )
      
      write.csv(
        removed_summary,
        file.path(output_dir, "rm_results", paste0("decontam_", group_name, ".csv")),
        row.names = FALSE
      )
    }
  }
  
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
  combined[, colSums(combined) > 0, drop = FALSE] |> as.matrix()
}

# ============================================================
# 2. Dataset-specific configuration
# ============================================================

decontam_groups <- list(
  Prokaryote = list(
    root1 = list(pattern = "Root", exclude = "G"),
    root2 = list(pattern = "G|pheno|ctab|freeze", exclude = NULL),
    soil1 = list(pattern = "Soil", exclude = "G"),
    soil2 = list(pattern = "G|ex", exclude = NULL)
  ),
  Fungi = list(
    root1 = list(pattern = "Root", exclude = "G"),
    root2 = list(pattern = "G_|pheno|ctab|freeze", exclude = NULL),
    root3 = list(pattern = "G2|pheno|ctab|freeze", exclude = NULL),
    root4 = list(pattern = "G3|pheno|ctab|freeze", exclude = NULL),
    root5 = list(pattern = "G4|pheno|ctab|freeze", exclude = NULL),
    soil  = list(pattern = "Soil", exclude = NULL)
  )
)

# ============================================================
# 3. Run
# ============================================================

for (dataset in names(decontam_groups)) {
  message("\n========================================")
  message("Decontamination: ", dataset)
  message("========================================")
  
  output_dir <- here("Output", dataset, "ASV_decontam")
  seqtab_nochim <- readRDS(file.path(output_dir, "merge_seqtab.nochim.rds"))
  
  dir.create(file.path(output_dir, "rm_results"), recursive = TRUE, showWarnings = FALSE)
  
  cleaned_tables <- list()
  
  for (group_name in names(decontam_groups[[dataset]])) {
    group <- decontam_groups[[dataset]][[group_name]]
    selected <- grepl(group$pattern, rownames(seqtab_nochim))
    
    group_seqtab <- seqtab_nochim[selected, , drop = FALSE]
    
    if (nrow(group_seqtab) == 0) {
      warning("No samples found for group: ", group_name)
      next
    }
    
    cleaned_tables[[group_name]] <- run_decontam(
      seqtab         = group_seqtab,
      neg_pattern    = "nega",
      remove_pattern = group$exclude,
      group_name     = group_name,
      output_dir     = output_dir
    )
  }
  
  seqtab_decontam <- combine_sequence_tables(cleaned_tables)
  saveRDS(seqtab_decontam, file.path(output_dir, "seqtab_decontam.rds"))
}
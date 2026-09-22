# 05_02_Alpha_diversity.R
#
# Purpose:
#   Calculate alpha diversity (observed richness) for microbial communities 
#   and perform statistical comparisons across habitats (Wilcoxon rank-sum test 
#   for Root/Soil) or sample types (GLMM with negative binomial distribution 
#   for Root&Soil).
#
# Input:
#   Output/03_Data_filtering/Seqdata/<data_type>/<sample_type>/
#     - <data_type>_<sample_type>_coverage_rared_th1fil.rds
#   Output/03_Data_filtering/Metadata/<data_type>/<sample_type>/
#     - <data_type>_<sample_type>_metadata_th1fil.csv
#
# Output:
#   Output/05_Community_characterization/02_Alfa_diversity/<data_type>/<sample_type>/
#     - <data_type>_<sample_type>_richness_alpha_df.rds
#     - <data_type>_<sample_type>_richness_stat_result.csv
#     - <data_type>_<sample_type>_richness_GLMM.rds (If Root&Soil)
#     - <data_type>_<sample_type>_richness_GLMM_summary.rds (If Root&Soil)
#     - <data_type>_<sample_type>_sample_size.csv (If Root&Soil)
#
# Analysis & Plot:
#   - Calculate species richness using vegan::specnumber.
#   - Conduct site-wise Wilcoxon rank-sum tests with BH p-value adjustment (Root or Soil).
#   - Fit negative-binomial GLMM (glmmTMB) accounting for random effects (Root&Soil).
#   - Perform model diagnostics via DHARMa and pairwise comparisons via emmeans.
#
# R version:
#   R 4.5.3
#
# Packages:
#   glmmTMB
#   DHARMa
#   emmeans
#   here
#   vegan
#   dplyr
#   rstatix
#   tidyr

# ======================================================================
# 1. Setup
# ======================================================================
library(glmmTMB)
library(DHARMa)
library(emmeans)
library(here)
library(vegan)
library(dplyr)
library(rstatix)
library(tidyr)

# Analysis settings -----------------------------------------------------
data_type <- "Fungi"   # "Prokaryote" or "Fungi"
sample_type <- "Root&Soil"  # "Root", "Soil", or "Root&Soil"

# Input and output directories ------------------------------------------
inputs <- here("Output", "03_Data_filtering", "Seqdata", data_type, sample_type)
inputm <- here("Output", "03_Data_filtering", "Metadata", data_type, sample_type)

output <- here("Output", "05_Community_characterization", "02_Alfa_diversity", data_type, sample_type)
dir.create(output, showWarnings = FALSE, recursive = TRUE)

# ======================================================================
# 2. Load data
# ======================================================================
seqdata <- readRDS(here(inputs, paste0(data_type, "_", sample_type, "_coverage_rared_th1fil.rds")))
metadata <- read.csv(here(inputm, paste0(data_type, "_", sample_type, "_metadata_th1fil.csv")))
rownames(metadata)<-metadata$Sample_ID

# ======================================================================
# 3. Calculate alpha diversity
# ======================================================================
alpha_df <- data.frame(
  Sample_ID = rownames(seqdata),
  Richness = specnumber(seqdata)
) |>
  left_join(
    metadata |>
      distinct(Sample_ID, site, habitat, sample_type, SID),
    by = "Sample_ID"
  ) |>
  mutate(
    habitat = factor(habitat, levels = c("Solfatara field", "Forest edge")),
    sample_type = factor(sample_type, levels = c("Root", "Soil"))
  )

# Save alpha-diversity data
saveRDS(
  alpha_df,
  file.path(output, paste0(data_type, "_", sample_type, "_richness_alpha_df.rds"))
)

# ======================================================================
# 4. Statistical analysis
# ======================================================================

# ----------------------------------------------------------------------
# Root or Soil only: Habitat comparison within each site (Wilcoxon test)
# ----------------------------------------------------------------------
if (sample_type %in% c("Root", "Soil")) {
  
  # Statistical test
  stats_df <- alpha_df |>
    group_by(site) |>
    rstatix::wilcox_test(Richness ~ habitat) |>
    ungroup() |>
    adjust_pvalue(method = "BH") |>
    mutate(
      group1 = "Solfatara field",
      group2 = "Forest edge",
      p.signif = case_when(
        p.adj < 0.001 ~ "***",
        p.adj < 0.01  ~ "**",
        p.adj < 0.05  ~ "*",
        TRUE          ~ "n.s."
      )
    )
  
  # Sample size
  n_df <- alpha_df |>
    filter(!is.na(Richness)) |>
    count(site, habitat, name = "n") |>
    pivot_wider(
      names_from = habitat,
      values_from = n,
      names_prefix = "n_"
    )
  
  # Effect size: Rank-biserial correlation
  effect_df <- alpha_df |>
    group_by(site) |>
    group_modify(~ {
      x <- .x$Richness[.x$habitat == "Solfatara field" & !is.na(.x$Richness)]
      y <- .x$Richness[.x$habitat == "Forest edge" & !is.na(.x$Richness)]
      
      if (length(x) < 1 || length(y) < 1) {
        return(tibble(effsize = NA_real_))
      }
      
      U <- wilcox.test(x, y, exact = FALSE)$statistic
      n1 <- length(x)
      n2 <- length(y)
      
      tibble(effsize = 2 * as.numeric(U) / (n1 * n2) - 1)
    }) |>
    ungroup()
  
  # Y position for plotting
  y_pos <- alpha_df |>
    group_by(site) |>
    summarise(
      y.position = max(Richness, na.rm = TRUE) * 1.05,
      .groups = "drop"
    )
  
  # Combine results
  stats_df <- stats_df |>
    left_join(n_df, by = "site") |>
    left_join(effect_df, by = "site") |>
    left_join(y_pos, by = "site")
  
  # ----------------------------------------------------------------------
  # Root & Soil: GLMM accounting for sampling-position dependence
  # ----------------------------------------------------------------------
} else if (sample_type == "Root&Soil") {
  
  # Negative-binomial GLMM
  model_nb <- glmmTMB(
    Richness ~ sample_type + (1 | SID),
    family = nbinom1,
    data = alpha_df
  )
  
  # Model diagnostics
  simulationOutput <- DHARMa::simulateResiduals(fittedModel = model_nb)
  dispersion_test <- DHARMa::testDispersion(simulationOutput)
  print(dispersion_test)
  
  # Model summary
  model_summary <- summary(model_nb)
  
  # Estimated marginal means
  emm <- emmeans(model_nb, ~ sample_type, type = "response")
  
  # Pairwise comparison (Root vs Soil)
  pairwise_df <- pairs(emm, reverse = TRUE) |>
    summary(infer = TRUE, type = "response") |>
    as.data.frame()
  
  # Sample size
  n_df <- alpha_df |>
    count(sample_type, name = "n") |>
    pivot_wider(
      names_from = sample_type,
      values_from = n,
      names_prefix = "n_"
    )
  
  # Number of sampling positions
  position_n_df <- alpha_df |>
    distinct(SID, sample_type) |>
    count(sample_type, name = "n_position") |>
    pivot_wider(
      names_from = sample_type,
      values_from = n_position,
      names_prefix = "n_position_"
    )
  
  # Significance annotation
  stats_df <- pairwise_df |>
    mutate(
      group1 = "Root",
      group2 = "Soil",
      p.signif = case_when(
        p.value < 0.001 ~ "***",
        p.value < 0.01  ~ "**",
        p.value < 0.05  ~ "*",
        TRUE            ~ "n.s."
      )
    ) |>
    left_join(n_df, by = character()) |>
    left_join(position_n_df, by = character())
  
  # Y position for plotting
  y_pos <- alpha_df |>
    summarise(y.position = max(Richness, na.rm = TRUE) * 1.05)
  
  stats_df <- stats_df |>
    mutate(y.position = y_pos$y.position)
  
} else {
  stop("Unexpected sample_type: ", sample_type, ". Use 'Root', 'Soil', or 'Root&Soil'.")
}

# Arrange common visualization columns
stats_df <- stats_df |>
  relocate(group1, group2, p.signif, y.position)

# ======================================================================
# 5. Save results
# ======================================================================
print(stats_df)

# Save statistical summary
write.csv(
  stats_df,
  file.path(output, paste0(data_type, "_", sample_type, "_richness_stat_result.csv")),
  row.names = FALSE
)

# Save GLMM objects and sample counts (Root&Soil only)
if (sample_type == "Root&Soil") {
  saveRDS(
    model_nb,
    file.path(output, paste0(data_type, "_", sample_type, "_richness_GLMM.rds"))
  )
  
  saveRDS(
    model_summary,
    file.path(output, paste0(data_type, "_", sample_type, "_richness_GLMM_summary.rds"))
  )
  
  metacount <- metadata |>
    group_by(site, habitat, sample_type) |>
    count()
  
  write.csv(
    metacount,
    file.path(output, paste0(data_type, "_", sample_type, "_sample_size.csv")),
    row.names = FALSE
  )
}
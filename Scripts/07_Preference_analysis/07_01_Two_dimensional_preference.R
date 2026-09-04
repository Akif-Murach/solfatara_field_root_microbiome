# 07_01_Two_dimensional_preference.R
#
# Purpose:
#   Perform Two-Dimensional Preference (2DP) analysis by randomizing sample assignments
#   and calculating observed vs. null distribution Z-scores and FDR-adjusted P-values.
#
# Input:
#   Loaded via "Scripts/07_Preference_analysis/07_00_Setup.R":
#     - seqdata: OTU / ASV abundance matrix
#     - metadata: Sample metadata
#     - focus: Target factor for preference analysis ("habitat", "host", etc.)
#     - nonfocus: Covariate factor to preserve during block sampling
#     - output: Directory path for saving results
#   Loaded via "Function/":
#     - blockSample.R: Function for block permutation preserving site structure
#     - Taxa.mat.R: Function to compute observed / randomized preference matrix
#
# Output:
#   Output directory (output):
#     - 2DP_Zvalue_<focus>.rds
#     - 2DP_two_sided_FDR_<focus>.rds
#
# Analysis:
#   - Parallelized block permutation using doParallel and doRNG.
#   - Z-score calculation based on unbiased variance of null distribution.
#   - Two-sided empirical P-value calculation and Benjamini-Hochberg (BH) FDR adjustment.
#
# R version:
#   R 4.5.3
#
# Packages:
#   here
#   foreach
#   doParallel
#   doRNG
#   stats

# ======================================================================
# 1. Setup
# ======================================================================
library(here)
# Load analysis settings -----------------------------------------------
source(here("Scripts", "07_Preference_analysis", "07_00_Setup.R"))
# Load functions -------------------------------------------------------
source(here("Function", "Blocksample.R"))
source(here("Function", "Taxa.mat.R"))

# Set random seed and parallel processing parameters -------------------
set.seed(1234)
n.core <- 9
rand <- 99999

# ======================================================================
# 2. Observed preference matrix
# ======================================================================
or_f2 <- Taxa.mat(
  t(seqdata),
  metadata,
  focus
)

nr <- nrow(or_f2)
nc <- ncol(or_f2)

# ======================================================================
# 3. Randomization (Parallel Processing)
# ======================================================================
# Initialize parallel backend ------------------------------------------
cl <- makeCluster(n.core)
registerDoParallel(cl)

clusterExport(
  cl,
  c(
    "seqdata",
    "metadata",
    "blockSample",
    "Taxa.mat",
    "or_f2"
  )
)

# Initialize result accumulators ---------------------------------------
init <- list(
  sum = matrix(0, nr, nc),
  sq_sum = matrix(0, nr, nc),
  pos_cnt = matrix(0L, nr, nc),
  neg_cnt = matrix(0L, nr, nc)
)

# Run Monte Carlo permutations -----------------------------------------
results <- foreach(
  r = 1:rand,
  .combine = function(a, b) {
    list(
      sum = a$sum + b$sum,
      sq_sum = a$sq_sum + b$sq_sum,
      pos_cnt = a$pos_cnt + b$pos_cnt,
      neg_cnt = a$neg_cnt + b$neg_cnt
    )
  },
  .init = init,
  .options.RNG = 1234
) %dorng% {
  
  # Randomize sample assignment while preserving site and non-focus grouping structure
  roc_f <- blockSample(
    seqdata,
    metadata[rownames(seqdata), "site"],
    metadata[rownames(seqdata), nonfocus],
    rownames(seqdata)
  )
  
  # Calculate preference matrix for the randomized data
  res_mat <- Taxa.mat(
    t(roc_f$matrix),
    metadata,
    focus
  )
  
  # Match rows and columns to the observed matrix
  res_mat <- res_mat[
    match(rownames(or_f2), rownames(res_mat)),
    match(colnames(or_f2), colnames(res_mat))
  ]
  
  stopifnot(
    identical(rownames(res_mat), rownames(or_f2)),
    identical(colnames(res_mat), colnames(or_f2))
  )
  
  list(
    sum = res_mat,
    sq_sum = res_mat^2,
    pos_cnt = (res_mat >= or_f2) & !is.na(res_mat),
    neg_cnt = (res_mat <= or_f2) & !is.na(res_mat)
  )
}

stopCluster(cl)

# ======================================================================
# 4. Z-scores calculation
# ======================================================================
rand_mean <- results$sum / rand

rand_var <- (results$sq_sum / rand) - rand_mean^2

# Convert to unbiased variance (finite sample size correction) ----------
rand_var <- rand_var * rand / (rand - 1)

rand_sd <- sqrt(rand_var)
rand_sd[rand_sd == 0] <- NA

z_mat <- (or_f2 - rand_mean) / rand_sd
z_mat[is.nan(z_mat)] <- NA

z_mat2 <- na.omit(z_mat)

# ======================================================================
# 5. Extract Z-scores and calculate two-sided P-values
# ======================================================================
# Because the two habitat categories are complementary, preference scores 
# and their corresponding Z-scores are perfectly anti-correlated; 
# therefore, only one habitat is retained for downstream analyses.

if (focus == "habitat") {
  
  z_A <- z_mat2[, "Solfatara field"]
  
  pos_cnt_A <- results$pos_cnt[, "Solfatara field"]
  neg_cnt_A <- results$neg_cnt[, "Solfatara field"]
  
  # Calculate two-sided P-values by doubling the smaller tail
  p_two <- 2 * pmin(
    (pos_cnt_A + 1) / (rand + 1),
    (neg_cnt_A + 1) / (rand + 1)
  )
  
  p_two[p_two > 1] <- 1
  
  p_two <- matrix(
    p_two,
    ncol = 1,
    dimnames = list(
      rownames(results$pos_cnt),
      "Solfatara field"
    )
  )
  
} else {
  
  # Calculate two-sided P-values by doubling the smaller tail
  p_upper <- (results$pos_cnt + 1) / (rand + 1)
  p_lower <- (results$neg_cnt + 1) / (rand + 1)
  
  p_two <- 2 * pmin(p_upper, p_lower)
  
  p_two[p_two > 1] <- 1
}

# ======================================================================
# 6. Multiple-testing correction (FDR)
# ======================================================================
# P-values are adjusted across OTUs using the Benjamini-Hochberg (BH) procedure.

p_vec <- as.vector(p_two)
valid <- !is.na(p_vec)

p_adj <- rep(NA_real_, length(p_vec))
p_adj[valid] <- p.adjust(p_vec[valid], method = "BH")

p_mat_fdr2 <- matrix(
  p_adj,
  nrow = nrow(p_two),
  ncol = ncol(p_two),
  dimnames = dimnames(p_two)
)

# ======================================================================
# 7. Save results
# ======================================================================
if (focus == "habitat") {
  
  saveRDS(
    z_A,
    file.path(
      output,
      paste0("2DP_Zvalue_", focus, ".rds")
    )
  )
  
} else {
  
  saveRDS(
    z_mat2,
    file.path(
      output,
      paste0("2DP_Zvalue_", focus, ".rds")
    )
  )
}

saveRDS(
  p_mat_fdr2,
  file.path(
    output,
    paste0("2DP_two_sided_FDR_", focus, ".rds")
  )
)
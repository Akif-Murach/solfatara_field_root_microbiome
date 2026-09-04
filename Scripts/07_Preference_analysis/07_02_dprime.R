# 07_02_dprime.R
#
# Purpose:
#   Calculate observed and null distribution of d' (d-prime) specialization index 
#   using block permutation tests. Computes Z-scores and FDR-adjusted P-values 
#   for interaction specialization (microbe or host level).
#
# Input:
#   Loaded via "Scripts/07_Preference_analysis/07_00_Setup.R":
#     - seqdata: OTU / ASV abundance matrix
#     - metadata: Sample metadata
#     - direction: Analysis direction ("microbe" or "host")
#     - output: Directory path for saving results
#   Loaded via "Function/":
#     - Blocksample.R: Function for block permutation preserving site and habitat structure
#     - Taxa.mat.R: Function to compute preference matrix
#     - (dfun): Function to compute d' specialization indices
#
# Output:
#   Output directory (output):
#     - dprime_Zvalue_<direction>.rds
#     - dprime_two_sided_FDR_<direction>.rds
#
# Analysis:
#   - Parallelized block permutation using doParallel and doRNG.
#   - Null distribution generation for d' values (permutations: 99,999 for microbe; 9,999 for host).
#   - Z-score computation and Benjamini-Hochberg (BH) FDR correction.
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
# Load analysis settings -----------------------------------------------
source(here("Scripts", "07_Preference_analysis", "07_00_Setup.R"))

# Load functions -------------------------------------------------------
source(here("Function", "Blocksample.R"))
source(here("Function", "Taxa.mat.R"))

# Set random seed and parallel processing parameters -------------------
set.seed(1234)
n.core <- 9

# ======================================================================
# 2. Observed preference matrix & d' calculation
# ======================================================================
or_f2 <- Taxa.mat(t(seqdata), metadata, "host")

if (direction == "microbe") {
  rand <- 99999
  d_prime_observed <- dfun(or_f2)
} else {
  rand <- 9999 # Reduced permutation iterations to optimize calculation cost
  d_prime_observed <- dfun(t(or_f2))
}

obs_d <- d_prime_observed$dprime
n_id <- length(obs_d)

# ======================================================================
# 3. Randomization (Parallel Processing)
# ======================================================================
# Initialize parallel backend ------------------------------------------
cl <- makeCluster(n.core)
registerDoParallel(cl)

mat_rows <- rownames(seqdata)
site_vec <- metadata[mat_rows, "site"]
habitat_vec <- metadata[mat_rows, "habitat"]

clusterExport(
  cl,
  c(
    "seqdata",
    "site_vec",
    "habitat_vec",
    "mat_rows",
    "metadata",
    "blockSample",
    "Taxa.mat",
    "dfun",
    "obs_d",
    "n_id"
  )
)

# Initialize result accumulators ---------------------------------------
init <- list(
  sum = numeric(n_id),
  sq_sum = numeric(n_id),
  pos_cnt = integer(n_id),
  neg_cnt = integer(n_id)
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
  
  # Randomize sample assignment ----------------------------------------
  random_matrix <- blockSample(
    seqdata,
    site_vec,
    habitat_vec,
    mat_rows
  )$matrix
  
  taxa_matrix <- Taxa.mat(t(random_matrix), metadata, "host")
  
  if (direction == "microbe") {
    d_rand <- dfun(taxa_matrix)$dprime
  } else {
    d_rand <- dfun(t(taxa_matrix))$dprime
  }
  
  # Ensure exact species ordering across iterations --------------------
  d_rand <- d_rand[names(obs_d)]
  stopifnot(identical(names(d_rand), names(obs_d)))
  
  list(
    sum = d_rand,
    sq_sum = d_rand^2,
    pos_cnt = (d_rand >= obs_d),
    neg_cnt = (d_rand <= obs_d)
  )
}

stopCluster(cl)

# ======================================================================
# 4. Z-scores & Empirical P-values calculation
# ======================================================================
# Mean and variance of null distribution -------------------------------
d_mean <- results$sum / rand

d_var <- (results$sq_sum / rand) - d_mean^2
d_var <- d_var * (rand / (rand - 1)) # Unbiased variance correction
d_sd <- sqrt(d_var)
d_sd[d_sd == 0] <- NA

# Z-scores -------------------------------------------------------------
z_d_prime <- (obs_d - d_mean) / d_sd

# Two-sided P-values (doubling the smaller tail) -----------------------
p_upper <- (results$pos_cnt + 1) / (rand + 1)
p_lower <- (results$neg_cnt + 1) / (rand + 1)
p_two <- 2 * pmin(p_upper, p_lower)
p_two[p_two > 1] <- 1

# ======================================================================
# 5. Multiple-testing correction (FDR)
# ======================================================================
p_value_fdr2 <- p.adjust(p_two, method = "BH")

# Inspect significant taxa ---------------------------------------------
which(p_value_fdr2 < 0.05)

# ======================================================================
# 6. Save results
# ======================================================================
saveRDS(
  z_d_prime,
  file.path(
    output,
    paste0("dprime_Zvalue_", direction, ".rds")
  )
)

saveRDS(
  p_value_fdr2,
  file.path(
    output,
    paste0("dprime_two_sided_FDR_", direction, ".rds")
  )
)
# 08_05_Phylogenetic_signal.R
# Purpose: Test phylogenetic signal in habitat and host preference using Blomberg's K and Pagel's lambda.
# Workflow: 1. Load tree & trait data -> 2. Prepare traits -> 3. Match OTUs -> 4. Calculate signal -> 5. Save results
# Input: Output/08_Phylogenetic_analysis/Phylo_data/
# Output: Output/08_Phylogenetic_analysis/Phylo_signal/

#===========================================================
# Packages & Setup
#===========================================================
library(here)
library(tidyverse)
library(ape)
library(phytools)

input_dir <- here("Output", "08_Phylogenetic_analysis", "Phylo_data")
input_dir2 <- here("Output", "07_Preference_analysis", "Proc_data",
                   "Fungi", "th3")
output_dir <- here("Output", "08_Phylogenetic_analysis", "Phylo_signal")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

#===========================================================
# 1. Load data & Prepare traits
#===========================================================
tree <- ape::read.tree(file.path(input_dir, "Newick_Pairwise_0.6.nwk"))
trait_data <- readRDS(file.path(input_dir2, "preference_data.rds")) |>
  as.data.frame() |>
  dplyr::rename(OTU = OTU_ID)

#===========================================================
# 2. Match OTUs between tree and trait data
#===========================================================
tree_tips <- tibble(tip.label = tree$tip.label) |>
  mutate(OTU = str_extract(tip.label, "F_[0-9]+"))

trait_tree_data <- tree_tips |>
  left_join(trait_data, by = "OTU") |>
  filter(!is.na(habitat_preference), !is.na(host_preference))

tree_trait <- drop.tip(tree, setdiff(tree$tip.label, trait_tree_data$tip.label))

#===========================================================
# 3. Calculate phylogenetic signal
#===========================================================
traits_to_test <- c("habitat_preference", "host_preference")
set.seed(1234)

phylogenetic_signal <- map_dfr(traits_to_test, \(trait_name) {
  trait_vector <- setNames(trait_tree_data[[trait_name]], 
                           trait_tree_data$tip.label) |> na.omit()
  tree_sub <- drop.tip(tree_trait, setdiff(tree_trait$tip.label, names(trait_vector)))
  
  result_K <- phylosig(tree_sub, trait_vector, method = "K", test = TRUE, nsim = 9999)
  result_lambda <- phylosig(tree_sub, trait_vector, method = "lambda", test = TRUE)
  
  tibble(
    trait = trait_name,
    method = c("Blomberg_K", "Pagel_lambda"),
    statistic = c(result_K$K, result_lambda$lambda),
    p_value = c(result_K$P, result_lambda$P)
  )
})

#===========================================================
# 4. Save results
#===========================================================
write.csv(phylogenetic_signal, 
          file = file.path(output_dir, "phylogenetic_signal_results.csv"), 
          row.names = FALSE)
phylogenetic_signal

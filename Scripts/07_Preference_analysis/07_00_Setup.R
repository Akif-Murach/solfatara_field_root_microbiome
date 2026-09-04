library(here)
library(tidyverse)
library(vegan)
library(parallel)
library(doParallel)
library(foreach)
library(doRNG)
library(bipartite)
#===========================================================
# Analysis settings
#===========================================================
data_type <- "Prokaryote"  # "Prokaryote" or "Fungi"
sample_type<-"Root" # only "Root"
analysis <- "dprime"       # "2DP" or "dprime"
focus <-  "habitat"        # "habitat" or "host"
nonfocus <- "host"    # "habitat" or "host"
direction <- "host"     # "microbe" or "host"
threshold <-5 # 1 or 3 or 5

#===========================================================
# Input
#===========================================================
inputs <- here("Output", "03_Data_filtering", "Seqdata", data_type, sample_type)
inputm <- here("Output", "03_Data_filtering", "Metadata", data_type, sample_type)

seq_abund <- readRDS(here(inputs, 
                          paste0(data_type, "_", sample_type,
                                 "_coverage_rared_th",threshold,"fil.rds")))

# Convert abundance to presence/absence
seqdata <- ifelse(seq_abund > 0,1,0)
# Safety check
all(seqdata %in% c(0, 1))


metadata <- read.csv(here(inputm, 
                          paste0(data_type, "_", sample_type,
                                 "_metadata_th",threshold,"fil.csv")),
                     row.names = 1)

#===========================================================
# Output
#===========================================================
output <- switch(
  analysis,
  "2DP" = here(
    "Output",
    "07_Preference_analysis",
　　analysis,
    data_type,
    focus,
    paste0("th",threshold)),
  "dprime" = here(
    "Output",
    "07_Preference_analysis",
    analysis,
    data_type,
    direction,
    paste0("th",threshold)))

dir.create(output, showWarnings = FALSE, recursive = TRUE)

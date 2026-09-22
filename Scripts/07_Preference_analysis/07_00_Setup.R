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
data_type <- "Fungi" # "Prokaryote" or "Fungi"
sample_type<-"Root" # only "Root"
# 2DP
focus <-  "habitat"        # "habitat" or "host"
nonfocus <- switch(
  focus,
  "habitat" = "host",
  "host" = "habitat",
  stop('focus must be "habitat" or "host".')
)
# dprime
direction <- "microbe"     # "microbe" or "host"

threshold <-3 # 1 or 3 or 5

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

# analysis and its output directory are set by 07_01 / 07_02.

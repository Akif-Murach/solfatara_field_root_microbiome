#06_00_Setup.R
#===========================================================
# Packages
#===========================================================
library(tidyverse)
library(vegan)
library(ape)
library(here)
library(ggnewscale)
library(ggstar)
#===========================================================
# Analysis settings
#===========================================================

data_type <- "Fungi" # "Prokaryote" or "Fungi"
sample_type <- "Root"     # "Root" or "Root&Soil"

#===========================================================
# Input
#===========================================================
inputs <- here("Output", "03_Data_filtering", "Seqdata", 
               data_type, sample_type)
inputm <- here("Output", "03_Data_filtering", "Metadata", 
               data_type, sample_type)

seqdata <- readRDS(
  here(inputs, paste0(
      data_type, "_", sample_type,
      "_coverage_rared_th3fil.rds")))

metadata <- read.csv(
  here(inputm, paste0(
      data_type, "_", sample_type,
      "_metadata_th3fil.csv")))
rownames(metadata)<-metadata$Sample_ID
#===========================================================
# Output
#===========================================================
output <- here("Output","06_Community_analysis",
               "00_Distance_matrix",data_type)
dir.create(output,showWarnings = FALSE,recursive = TRUE)

#===========================================================
# Distance matrix
#===========================================================
# Sørensen dissimilarity index
dist <- vegdist(seqdata, method = "bray", binary=TRUE)

saveRDS(dist,file=here(output, 
                       paste0(sample_type, "_", "distance_matrix.rds")))

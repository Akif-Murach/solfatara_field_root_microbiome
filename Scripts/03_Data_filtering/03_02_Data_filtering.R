# 03_02_Data_filtering.R
library(here)

#===========================================================
# Argument settings
#===========================================================

data_type <- "Prokaryote"   # "Prokaryote" or "Fungi"
sample_type <- "Root&Soil"       # "Root" or "Soil" or "Root&Soil"
th <- 3 # OTU occurence threshold >= 1 or 3 or 5 
#===========================================================
# Input file
#===========================================================
# Read rarefied OTU table and metadata

input1 <- here("Output", "01_Data_processing", "Covrfy", data_type)
input1rs<-here("Output", "03_Data_filtering", "Seqdata", data_type,
               "Root&Soil")
#OTU table
seqdata <- readRDS(here(
  ifelse(sample_type=="Root&Soil",input1rs,input1),
  paste0(data_type, "_", sample_type,"_coverage_rared.rds")))

#metadata 
input2<- here("Output", "02_Plant_root_identification","Metadata")

metadata <- read.csv(here(input2,data_type,
  paste0(data_type, "_", sample_type,"_metadata.csv")))
rownames(metadata)<-metadata$Sample_ID
#===========================================================
# Output directory
#===========================================================

output1 <- here("Output", "03_Data_filtering", "Seqdata",
              data_type, sample_type)
dir.create(output1,showWarnings = FALSE,recursive = TRUE)
output2<- here("Output", "03_Data_filtering", "Metadata",
               data_type, sample_type)
dir.create(output2,showWarnings = FALSE,recursive = TRUE)

#===========================================================
# Match samples and rare OTU filtering 
#===========================================================
common_samples <- intersect(rownames(metadata), rownames(seqdata))

seqdataf <- seqdata[common_samples, ]
metaf    <- metadata[common_samples, ]

# Filter OTUs appearing in at least 3 samples
otu_keep   <- colSums(seqdataf > 0) >= th
seqdata_th <- seqdataf[, otu_keep]

# Filter valid samples that retain at least one OTU
sample_keep <- rowSums(seqdata_th > 0) > 0
df_fil      <- seqdata_th[sample_keep, ]
meta_fil    <- metaf[sample_keep, ]

# Check matrix dimensions (samples x OTUs)
dim(df_fil)

#===========================================================
# Saving
#===========================================================
saveRDS(df_fil,
        file = here(
          output1,
          paste0(
            data_type, "_", sample_type,
            paste0("_coverage_rared_th",th,"fil.rds"))))

write.csv(meta_fil,
  file = here(
  output2,
  paste0(
    data_type, "_", sample_type,
    paste0("_metadata_th",th,"fil.csv"))))


library(here)

data_type <- "Fungi"   # "Prokaryote" or "Fungi"
input <- here("Output", "01_Data_processing", "Covrfy", data_type)

output<-here("Output", "03_Data_filtering", "Seqdata",
             data_type, "Root&Soil")
dir.create(output,showWarnings = FALSE,recursive = TRUE)


root_seq <- readRDS(here(
  input,paste0(data_type,"_Root_coverage_rared.rds")))
soil_seq <- readRDS(here(
  input,paste0(data_type,"_Soil_coverage_rared.rds")))

mergedata <- bind_rows(as.data.frame(root_seq), 
                       as.data.frame(soil_seq))

mergedata[is.na(mergedata)] <- 0

saveRDS(mergedata, 
        file.path(output, 
                  paste0(data_type,"_Root&Soil_coverage_rared.rds")))

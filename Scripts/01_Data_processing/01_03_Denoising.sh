#!/bin/bash

####################################################################
## 											
## ---------- 	   Denoising sequence by DADA2      ------------- ##
##
####################################################################

inputdir=03_FilterTrimming_fastaFiles
outputdir=04_Denoising
thread=8 #- Please set the number of threads to use
minident=1
vsearchpath=($ which vsearch) #- Please set the path of VSEARCH
## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ##
piplineID="04_Denoising"

####################################################################

## ============= Remove and Make directories ==================== ##

## --Remove older directory

if [ -d ${piplineID} ]; then
	rm -r ${piplineID}
fi

if ls log_and_script/${piplineID}_log* >/dev/null 2>&1 ; then
	rm log_and_script/log${piplineID}
fi

## ------------------------------------------------------------- ##
## -- Making directory to save results
mkdir -p ${piplineID}

## ------------------------------------------------------------- ##
## -- Version check
echo "## ++++++ ${piplineID} +++++ ##" >> log_and_script/Version.txt

####################################################################

cat <<RRR > log_and_script/script${piplineID}.R

########################################################################
## -- DADA2
## -- Running on R

## == Setting parameter
inputdir <- "$inputdir" 
outputdir <- "$outputdir"
minident=${minident}
vsearchpath="${vsearchpath}"

## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ##
piplineID="${piplineID}"

start <- Sys.time()
print(start); cat("\n")
########################################################################

RRR

cat <<'RRR' >> log_and_script/script${piplineID}.R

library(dada2)
library(seqinr)

dir.create(outputdir, showWarnings = FALSE)

takeSamplenames <- function(x, sep=''){
  split <- do.call(rbind, strsplit(x, sep) )
  check.unique <- sapply(apply(split, 2, table), length)
  uniquename <- which(check.unique>1)
  if(length(uniquename)==1  ) {
    samplename <- split[,uniquename]
  }else{
    samplename <- apply(split[,uniquename], 1, paste, collapse="_")
  }
  return(samplename)
}

########################################################################

run_dirs <- list.dirs(inputdir, recursive = FALSE)

if (length(run_dirs) == 0) {
  stop("The path directory is missing\n")
}

########################################################################

for (path in run_dirs) {

  run_name <- basename(path)
  cat("Processing:", run_name, "\n")

  fnFs <- sort(list.files(path, pattern="fastq.gz", full.names = TRUE))

  if (length(fnFs) == 0) {
    cat("No files in", run_name, "\n")
    next
  }

  # Output for each sequencing run
  run_outdir <- file.path(outputdir, run_name)
  dir.create(run_outdir, showWarnings = FALSE, recursive = TRUE)

  cat( sprintf("Start learnerror from %s...\n", Sys.time()))

  errF <- invisible( learnErrors(fnFs, nbases=1e8, multithread=TRUE) )

  pdf(sprintf('%s/plotErrors.pdf', run_outdir))
  print(plotErrors(errF, nominalQ=TRUE))
  dev.off()

  derepFs <- derepFastq(fnFs, verbose = TRUE)

  samplenames <- takeSamplenames(basename(fnFs), sep='__')

  names(derepFs) <- samplenames

  dadaFs <- invisible( dada(derepFs, err = errF, multithread = TRUE) )

  st.all <- makeSequenceTable(dadaFs)

  saveRDS(st.all, sprintf('%s/%s_stall_no_rmchimera.rds', run_outdir, run_name))

  cat("Finished:", run_name, "\n\n")
}

##################################################################################################

finish <- Sys.time()
print(finish-start); cat("\n")

##################################################################################################
sink("log_and_script/Version.txt", append=TRUE)
sessionInfo()
sink()

RRR

Rscript log_and_script/script${piplineID}.R 2>> log_and_script/log${piplineID}.txt
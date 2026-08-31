#!/bin/bash

####################################################################
## ---------- DADA2 filtering only ---------- ##
####################################################################
## Input directory
fastadir=02_Cutadaptor_fastaq

## -- Options
minLen=200
minQ=10
trancWindow=5
maxEE=3
truncR=240
avgqual=15

thread=8 #- Number of threads

piplineID="03_FilterTrimming"

## -- Version check (optional)
echo "## ++++++ ${piplineID} +++++ ##" >> log_and_script/Version.txt

####################################################################
echo "Start at $timestamp"

cat <<RRR > log_and_script/script${piplineID}.R

########################################################################
## -- DADA2
## -- Running on R

## == Setting parameter
minLen=$minLen
minQ=$minQ
maxEE=$maxEE
truncR=$truncR

inputdir="$fastadir"

piplineID="$piplineID"

########################################################################

RRR

cat <<'RRR' >> log_and_script/script${piplineID}.R

# -- Load library
library(dada2) 

run_dirs <- list.dirs(inputdir, recursive = FALSE)

if (length(run_dirs) == 0) {
  stop("The path directory is missing\n")
}

## -- Set random seed
ran.seed <- 1234
set.seed(ran.seed)

## -- Filtering for each run
for (path in run_dirs) {

  fnFs <- sort(list.files(path, pattern="fastq.gz", full.names = TRUE))

  if (length(fnFs) == 0) next

  ## -- Output directory
  outdir <- file.path(sprintf("%s_fastaFiles", piplineID), basename(path))
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

  filtFs <- file.path(outdir, basename(fnFs))

  ## -- Filtering
  filterAndTrim(fnFs, filtFs, maxN = 0, 
                         maxEE = maxEE, minLen = minLen, truncLen=truncR, minQ=minQ, truncQ=11,
                         rm.phix = TRUE, compress = TRUE, multithread = TRUE) 
}

sink("log_and_script/Version.txt", append=TRUE)
sessionInfo()
sink()

RRR

Rscript log_and_script/script${piplineID}.R

## -- Remove tiny files
find ${piplineID}_fastaFiles -size -500c -delete

echo "Finish at $timestamp"
####################################################################
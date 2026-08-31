#!/bin/bash

####################################################################
## 											
## ---------- 	  Check and cut sequence adaptor     ------------ ##
## 
####################################################################

inputdir=01_Demultiplexed_fastaq #- Please input download data from DDBJ
outputdir=02_Cutadaptor_fastaq
fseq=NNNNNNGTGYCAGCMGCCGCGGTAA #- 16S 515f Prokaryotes 
rseq=NNNNNNGGACTACNVGGGTWTCTAAT #- 16S 806rB Prokaryotes

fseq=NNNNNNCTHGGTCATTTAGAGGAASTAA　#- ITS1F_KYO1 Fungi 
rseq=NNNNNNTTYRCTRCGTTCTTCATC　#- ITS2_KYO2 Fungi

fseq=NNNNNNACBTRGTGTGAATTGCAGRATC　#- ITS_3p62plF1 Plant
rseq=NNNNNNTCCTCCGCTTATTKATATGC　#- ITS_4unR1 Plant

cutadaptpath=$(which cutadapt) #- Please set the path of Cutadapt
## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ##
piplineID=02_Cutadaptor

####################################################################

## ============= Remove and Make directories ==================== ##

## -- Remove older directory
if [ -d $outputdir ]; then
  rm -r $outputdir
fi

## -- Version check
echo "## ++++++ ${piplineID} +++++ ##" >> log_and_script/Version.txt
cv=`cutadapt --version`
echo "cutadapt ${cv}" >> log_and_script/Version.txt

########################################################################
## -- Cutadaptor by R

cat <<RRR > log_and_script/script${piplineID}.R

inputdir="${inputdir}"
outputdir="${outputdir}"
primerf="${fseq}"
primerr="${rseq}"
thread=${thread}
cutadaptpath="${cutadaptpath}"

## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ##

piplineID="${piplineID}"

start <- Sys.time()
print(start); cat("\n")
########################################################################
RRR

cat <<'RRR' >> log_and_script/script02_Cutadaptor.R
## ===================== Definition of function =====================  ##
primerHits <- function(primer, fn) {
   # Counts number of reads in which the primer is found
   nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
    return(sum(nhits > 0))
  }

allOrients <- function(primer) {
    # Create all orientations of the input sequence
    require(Biostrings)
    dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
    orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
    return(sapply(orients, toString))  # Convert back to character vector
}

Check.primer = function(fwd= fwd, rev,
                        path=path, input.fastq=fnFs){

  FWD.orients <- allOrients(fwd)
  REV.orients <- allOrients(rev)

  files <- list.files(path, full.names = TRUE)

  tmp1 <- sapply(FWD.orients, function(p) {
     sum(sapply(files, primerHits, primer = p))
  })
  tmp2 = sapply(REV.orients, primerHits, fn = list.files(path, full.names = TRUE))
  rbind(tmp1,  tmp2)
  
}

########################################################################
######  ====================== Main part ======================  #######
########################################################################

# -- Create base output directory
dir.create(outputdir, showWarnings = FALSE)

# -- Load library
library(seqinr); library(stringr); library(ShortRead)
library(Biostrings); library(dada2); library(doParallel)

run_dirs <- list.dirs(inputdir, recursive = FALSE)

if (length(run_dirs) == 0) {
  stop("No Run directories found")
}

for (run_path in run_dirs) {

  run_name <- basename(run_path)
  cat("Processing:", run_name, "\n")

  # Output folder for each sequencing run
  run_outdir <- file.path(outputdir, run_name)
  dir.create(run_outdir, showWarnings = FALSE)

  fnFs <- sort(list.files(run_path, pattern="\\.fastq.gz$", full.names = TRUE))

  if (length(fnFs) == 0) {
    cat("No FASTQ files in", run_name, "\n")
    next
  }

  # Checking primers
  FWD.ForwardReads = Check.primer(
    fwd = primerf,
    rev = primerr,
    input.fastq = fnFs,
    path = run_path
  )

  cat(sprintf('Primer hits in %s: %s\n', run_name, sum(FWD.ForwardReads)))

  print(FWD.ForwardReads)

  # Output file name
  fnFs.cut <- file.path(run_outdir, basename(fnFs))

  RVS.RC <- dada2:::rc(primerr)
  R1.flags <- paste("-g", primerf, "-a", RVS.RC)

  cl <- makeCluster(detectCores(logical = FALSE))
  registerDoParallel(cl)

  tmp <- foreach(f = seq_along(fnFs)) %dopar% {
    invisible(system2(
      cutadaptpath,
      args = c(
        R1.flags,
        "-n", 2,
        "-j", 2,
        "--max-n", 0,
        "-o", fnFs.cut[f],
        "-m", 10,
        fnFs[f]
      )
    ))
  }

  stopCluster(cl)

  cat("Finished:", run_name, "\n\n")
}

finish <- Sys.time()
print(sprintf("Finish. %s", finish-start)); cat("\n")

RRR
Rscript log_and_script/script02_Cutadaptor.R 2>&1 | tee  log_and_script/log02_Cutadaptor.txt

find 02_Cutadaptor_fastaq -size -500c -delete


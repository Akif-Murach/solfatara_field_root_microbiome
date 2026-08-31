#!/bin/bash
set -e

inputdir=Summary
outputdir=05_TaxonomyAnnotation
thread=32
referencepath=overall_genus
minident=0.97
maxpopposer=0.05
minsoratio=19
NNC=5,90%

mkdir -p ${outputdir}
mkdir -p log_and_script

clmakecachedb \
--blastdb=${referencepath} \
${inputdir}/merge_clus_seq.fasta \
${outputdir}/01_cachedb_species

clidentseq \
--method=QC \
--blastdb=${outputdir}/01_cachedb_species \
${inputdir}/merge_clus_seq.fasta \
${outputdir}/02_neighborhoods_QC.txt

classigntax \
--taxdb=${referencepath} \
${outputdir}/02_neighborhoods_QC.txt \
${outputdir}/03_taxalist_QC_strict.tsv

classigntax \
--taxdb=${referencepath} \
--maxpopposer=${maxpopposer} \
--minsoratio=${minsoratio} \
${outputdir}/02_neighborhoods_QC.txt \
${outputdir}/03_taxalist_QC_relax.tsv

clidentseq \
--method=${NNC} \
--blastdb=${outputdir}/01_cachedb_species \
${inputdir}/merge_clus_seq.fasta \
${outputdir}/02_neighborhoods_NNC.txt

classigntax \
--taxdb=${referencepath} \
--minnsupporter=1 \
${outputdir}/02_neighborhoods_NNC.txt \
${outputdir}/03_taxalist_NNC.tsv

clmergeassign \
--preferlower \
--priority=descend \
${outputdir}/03_taxalist_QC_strict.tsv \
${outputdir}/03_taxalist_QC_relax.tsv \
${outputdir}/03_taxalist_NNC.tsv \
${outputdir}/taxonomy_merged.tsv

clfillassign \
${outputdir}/taxonomy_merged.tsv \
${outputdir}/taxonomy_merged_filled.tsv

# =========================
# R reformat
# =========================
cat <<RRR > log_and_script/script05_claident_reformat.R
inputdir="\${inputdir}"
minident=\${minident}
outputdir="\${outputdir}"

library(seqinr)
library(stringr)
library(dplyr)
library(tidyr)

df <- read.table(sprintf("%s/taxonomy_merged_filled.tsv",outputdir), header=T, row.names=1, sep="\t")

otu_asv <- read.table(sprintf("%s/ASV_OTU_corestab_%s.txt", inputdir, minident),
header=TRUE, row.names=2)[,-c(1:2)]

otu_asv_lf = cbind(ASV=rownames(otu_asv), otu_asv) |>
  pivot_longer(-1, names_to="OTU") |>
  filter(value>0)

seqfasta <- read.fasta(sprintf("%s/merge_clus_seq.fasta",inputdir))
seqlist <- sapply(seqfasta, paste, collapse="")

df2 <- df[, intersect(colnames(df),
c("superkingdom","kingdom","phylum","class","order","family","genus","species"))]

colnames(df2) <- stringr::str_to_title(colnames(df2))

mat <- as.matrix(df2)
mat[mat==""] <- "Unidentified"
mat <- gsub(" ", "_", mat)

df3 <- as.data.frame(mat)
df3[is.na(df3)] = "Unidentified"

taxa_list = left_join(otu_asv_lf,
data.frame(ASV=rownames(df3), df3, seq=seqlist[rownames(df3)]),
by="ASV")

write.table(cbind(ID=rownames(df3), df3),
sprintf("%s/taxonomy_list.txt",outputdir),
row.names=F, quote=F, sep="\t")

saveRDS(cbind(ID=rownames(df3), df3),
sprintf("%s/taxonomy_list.rds",outputdir))

RRR

Rscript log_and_script/script05_claident_reformat.R

chmod +x log_and_script/script05_claident.sh

echo "Script ready:"
echo "bash log_and_script/script05_claident.sh"
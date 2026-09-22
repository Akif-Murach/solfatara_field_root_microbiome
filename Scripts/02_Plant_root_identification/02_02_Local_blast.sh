#!/usr/bin/env bash
set -euo pipefail

# 02_02_Local_blast.sh
#
# Identify the most similar reference-leaf sequences for Plant root OTUs.
# Run from the repository root after 02_01_Process_seq_and_taxa.R.

threads=8
evalue=1e-20

input_dir="Output/02_Plant_root_identification/Seqdata"
output_dir="Output/02_Plant_root_identification/Local_Blast"
leaf_fasta="${input_dir}/Leaf.fasta"
root_fasta="${input_dir}/Root.fasta"
db_prefix="${output_dir}/Leaf_db"
blast_output="${output_dir}/Leaf_Root_blast.tsv"

mkdir -p "${output_dir}"

# Keep every reference subject eligible for reporting, rather than only five.
max_target_seqs=$(awk '/^>/ {n++} END {print n+0}' "${leaf_fasta}")
if (( max_target_seqs < 1 )); then
  echo "No reference sequences found in ${leaf_fasta}" >&2
  exit 1
fi

makeblastdb \
  -in "${leaf_fasta}" \
  -dbtype nucl \
  -out "${db_prefix}"

blastn \
  -query "${root_fasta}" \
  -db "${db_prefix}" \
  -out "${blast_output}" \
  -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" \
  -max_target_seqs "${max_target_seqs}" \
  -evalue "${evalue}" \
  -num_threads "${threads}"

echo "BLAST output: ${blast_output}"

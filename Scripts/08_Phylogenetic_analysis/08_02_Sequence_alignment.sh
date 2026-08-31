#!/bin/bash

# ===========================================================
# 05_02_Sequence_alignment.sh
#
# Purpose:
#   Multiple sequence alignment and alignment trimming
#   for the Hyaloscyphaceae phylogenetic analysis.
#
# Software:
#   MAFFT v7.526
#   trimAl v1.5.1
# ===========================================================


#===========================================================
# Input / Output
#===========================================================

input_dir="Output/05_Phylogenetic_analysis/Phylo_data"


#===========================================================
# 1. Multiple sequence alignment
#===========================================================

# MAFFT v7.526
# E-INS-i strategy with 16 iterative refinements

mafft \
  --genafpair \
  --maxiterate 16 \
  --inputorder \
  "${input_dir}/merged_Hyaloscyphaceae_Ref.fasta" \
  > "${input_dir}/refseq_amplicon_merge_align.fasta"


#===========================================================
# 2. Convert sequence characters to uppercase
#===========================================================

# trimAl returned an error when lowercase "n" characters were present.
# Sequence headers are retained unchanged.

awk '
/^>/ {print; next}
{print toupper($0)}
' \
  "${input_dir}/refseq_amplicon_merge_align.fasta" \
  > "${input_dir}/refseq_amplicon_merge_align_upper.fasta"


#===========================================================
# 3. Alignment trimming
#===========================================================

# trimAl v1.5.1
#
# The final phylogenetic analysis used a gap threshold of 0.6.
# -gt 0.6: retain positions with at least 60% non-gap characters.

trimal \
  -in "${input_dir}/refseq_amplicon_merge_align_upper.fasta" \
  -out "${input_dir}/trimmed_0.6.fasta" \
  -gt 0.6
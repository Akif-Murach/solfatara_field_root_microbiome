#!/usr/bin/env bash
set -e

# 01_07_Claident_plant_annotation.sh
#
# Purpose:
#   Assign Plant taxonomy to every ASV with Claident and reformat the
#   assignments. OTU representatives are extracted in 01_08.
#
# Usage:
#   bash Scripts/01_Data_processing/01_07_Claident_plant_annotation.sh \
#     /path/to/claident/database_prefix
#
# Alternatively, set CLAIDENT_DB before running the script.

# ======================================================================
# 1. Paths and settings
# ======================================================================
script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
project_root=${PROJECT_ROOT:-$(cd -- "${script_dir}/../.." && pwd)}

input_dir="${project_root}/Output/01_Data_processing/Plant"
output_dir="${input_dir}/Claident"
reference_db=${1:-${CLAIDENT_DB:-}}

maxpopposer=0.05
minsoratio=19
nnc="5,90%"

if [ -z "${reference_db}" ]; then
  echo "Provide the Claident database prefix as the first argument or set CLAIDENT_DB." >&2
  exit 1
fi

input_fasta="${input_dir}/merge_clus_seq.fasta"
if [ ! -f "${input_fasta}" ]; then
  echo "Input FASTA not found: ${input_fasta}" >&2
  exit 1
fi

mkdir -p "${output_dir}"

# ======================================================================
# 2. Claident taxonomy assignment
# ======================================================================
clmakecachedb \
  --blastdb="${reference_db}" \
  "${input_fasta}" \
  "${output_dir}/01_cachedb_species"

clidentseq \
  --method=QC \
  --blastdb="${output_dir}/01_cachedb_species" \
  "${input_fasta}" \
  "${output_dir}/02_neighborhoods_QC.txt"

classigntax \
  --taxdb="${reference_db}" \
  "${output_dir}/02_neighborhoods_QC.txt" \
  "${output_dir}/03_taxalist_QC_strict.tsv"

classigntax \
  --taxdb="${reference_db}" \
  --maxpopposer="${maxpopposer}" \
  --minsoratio="${minsoratio}" \
  "${output_dir}/02_neighborhoods_QC.txt" \
  "${output_dir}/03_taxalist_QC_relax.tsv"

clidentseq \
  --method="${nnc}" \
  --blastdb="${output_dir}/01_cachedb_species" \
  "${input_fasta}" \
  "${output_dir}/02_neighborhoods_NNC.txt"

classigntax \
  --taxdb="${reference_db}" \
  --minnsupporter=1 \
  "${output_dir}/02_neighborhoods_NNC.txt" \
  "${output_dir}/03_taxalist_NNC.tsv"

clmergeassign \
  --preferlower \
  --priority=descend \
  "${output_dir}/03_taxalist_QC_strict.tsv" \
  "${output_dir}/03_taxalist_QC_relax.tsv" \
  "${output_dir}/03_taxalist_NNC.tsv" \
  "${output_dir}/taxonomy_merged.tsv"

clfillassign \
  "${output_dir}/taxonomy_merged.tsv" \
  "${output_dir}/taxonomy_merged_filled.tsv"

# ======================================================================
# 3. Reformat ASV taxonomy
# ======================================================================
Rscript --vanilla - "${input_dir}" "${output_dir}" <<'RSCRIPT'
args <- commandArgs(trailingOnly = TRUE)
input_dir <- args[[1]]
output_dir <- args[[2]]

taxonomy <- read.table(
  file.path(output_dir, "taxonomy_merged_filled.tsv"),
  header = TRUE,
  row.names = 1,
  sep = "\t",
  check.names = FALSE,
  quote = "",
  comment.char = ""
)

ranks <- c(
  "superkingdom", "kingdom", "phylum", "class",
  "order", "family", "genus", "species"
)
taxonomy <- taxonomy[, intersect(colnames(taxonomy), ranks), drop = FALSE]
colnames(taxonomy) <- tools::toTitleCase(colnames(taxonomy))

taxonomy <- as.matrix(taxonomy)
taxonomy[taxonomy == ""] <- "Unidentified"
taxonomy <- gsub(" ", "_", taxonomy)
taxonomy[is.na(taxonomy)] <- "Unidentified"

taxonomy <- as.data.frame(taxonomy, stringsAsFactors = FALSE)
taxonomy_table <- cbind(ID = rownames(taxonomy), taxonomy)

write.table(
  taxonomy_table,
  file.path(output_dir, "taxonomy_list.txt"),
  row.names = FALSE,
  quote = FALSE,
  sep = "\t"
)
saveRDS(taxonomy_table, file.path(output_dir, "taxonomy_list.rds"))
saveRDS(taxonomy, file.path(input_dir, "merge_taxonomylist.rds"))
RSCRIPT

echo "Plant Claident annotation completed: ${output_dir}"

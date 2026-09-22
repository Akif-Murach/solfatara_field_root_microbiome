# Fumarole-field microbiome analysis

Scripts for **Environmental filtering and host identity collectively shape
root-associated microbiomes of Ericaceae and ectomycorrhizal plants in
fumarole fields**.

## Workflow at a glance

```text
Raw reads (DDBJ: PRJDB42913)
  01_01–01_03  Primer removal, QC and denoising              [Linux; per run]
       ↓ transfer run-level sequence tables
  01_04–01_05  Merge runs, remove chimeras and contaminants  [Mac]
       ↓
  01_06–01_08  Taxonomy assignment and 97% OTU clustering    [Linux]
       ↓
  01_09        Normalize sample IDs                        [Mac]
       ↓
Processed OTU tables, taxonomy and representative sequences (Data/)
       ├── 01_10  Microbial rarefaction ──────────┐
       └── 02     Plant-root identification ─────┤
                                                ↓
                                      03  Community filtering
                                                ├── 05  Composition and alpha diversity
                                                ├── 06  Community-level analyses
                                                └── 07  Habitat/host preference
                                                          ↓
                                                    08  Fungal phylogeny
Soil chemistry data ─────────────────────── 04  Soil analyses
```

Data pre-processing steps (01_01–01_03 and 01_06–01_08) were run on 
`x86_64-pc-linux-gnu` under `Ubuntu 20.04.6 LTS`.
All steps other than the above were run on 
`aarch64-apple-darwin20` under `macOS Tahoe 26.6.2`.
The sections below describe the two entry points (Route (A) and (B)) 
and the necessary handoffs;
individual scripts contain the detailed settings and methods.

# Route (A): ###################################################################
## Start from processed data

Open `Solfatara_Analysis.Rproj` and retain the supplied `Data/` hierarchy.
For each of `Prokaryote`, `Fungi` and `Plant`, the principal inputs in
`Data/<dataset>/Seqdata/` are `seqOTUtab_filtered.rds`,
`OTU_merge_taxonomylist.rds` and `OTUseq_0.97.fasta`.
Retain the accompanying metadata, `Data/Soil_analysis/`,
`Data/Plant/Supplementary_annotation/LBlast_noref_ano.csv` and
`Data/Fungi/Seqdata/Hyaloscyphaceae_refseq.fasta`.
Tables, taxonomy, sequences and curated annotations must belong to the same OTU release.

1. **Skip 01_01–01_09.** Run `01_10_Coverage_rarefaction.R` for each of
   Prokaryote/Fungi × Root/Soil by editing `data_type` and `sample_type`.
   Plant data do not pass through 01_10.
2. Run stage **02** in numerical order to generate host metadata. This stage
   requires BLAST+ even when starting from processed data. For a new OTU
   reconstruction, curate the unmatched OTUs in `LBlast_noref.csv`, prepare
   the corresponding `LBlast_noref_ano.csv`, and rerun 02_03 before 02_04.
3. Run **03_01** once per microbial marker to merge Root and Soil data, then
   **03_02** for each required dataset, sample type and occurrence threshold.
   Use th1 for alpha diversity, th3 for the main community/preference analyses,
   and th1/th3/th5 for preference robustness comparisons.
4. Continue with stages **04–08** as needed. Run analyses before their figure
   and supplementary-table scripts.

## Notes

- Stage 03 writes matching sequence tables and metadata below
  `Output/03_Data_filtering/{Seqdata,Metadata}/`. Regenerate both together when
  upstream inputs change.

- Settings are edited in the scripts, not supplied as command-line arguments
  unless explicitly supported. Stage 06 uses `06_00_Setup.R`; stage 07 uses
  `07_00_Setup.R`. Generate each dataset/threshold/focus/direction combination
  needed by the downstream figures and comparisons. The 07_01 and 07_02 entry
  points select 2DP and d-prime, respectively.

- Phylogenetic tree reconstruction
  After the required fungal preference results have been generated, run
  08_01 (sequence selection), then 08_02 (MAFFT/trimAl). 
  Reconstruct the tree in MEGA using [`08_03_MEGA12_NJ_tree_settings.txt`]
  (Scripts/08_Phylogenetic_analysis/08_03_MEGA12_NJ_tree_settings.txt).
  Save the tree as `Output/08_Phylogenetic_analysis/Phylo_data/Newick_Pairwise_0.6.nwk`
  before running 08_04 and 08_05.


# Route (B): ###################################################################
## Start from raw reads

Raw FASTQ files are available through DDBJ BioProject **PRJDB42913** and are
not stored here. Use its run records to obtain the demultiplexed reads.

### Run-level processing: 01_01–01_03 (Linux)

Create a separate working directory for each marker and sequencing run,
containing `01_Demultiplexed_fastaq/` and `log_and_script/`.
Place only the selected `.fastq.gz` files in the former. Preserve complete
sample IDs, including control and resequencing tokens (`nega`, `_G_`, `_G2_`).
Use a common filename suffix separated from the sample ID by `__` and check
the sample names in the resulting sequence table.

**01_01–01_03 generate scripts; they do not complete the processing alone.**
Run each generator with its documented arguments, then execute its generated
script in the same run directory:

| Generator | Generated script to execute |
|---|---|
| `01_01_cutadapt_R.sh` | `log_and_script/script02_Cutadaptor.sh` |
| `01_02_readQC_dada2.sh` | `log_and_script/script03_FilterTrimming.sh` |
| `01_03_Denoising_to_prechimera.sh` | `log_and_script/script04_Denoising.sh` |

Use the marker-specific values in
[`parameterList.tsv`](Scripts/01_Data_processing/parameterList.tsv);
the wrappers do not read this file automatically. Use a fresh working
directory for repeats because generated scripts replace derived outputs.
Check logs, retained samples and controls before proceeding.

Transfer each resulting `04_Denoising/stall_no_rmchimera.rds` to
`Data/<dataset>/Seqdata/<run>/stall_no_rmchimera.rds` on the Mac.
The required run-directory names are listed in `01_04_Merge_runs.R`.

### Combined processing: 01_04–01_09

Run from the repository root, carrying the required `Data/` and `Output/`
files between machines while preserving their relative paths:

1. **Mac — 01_04, 01_05:** merge runs, remove chimeras and perform
   decontamination. Preserve negative controls and resequencing identifiers
   until this step; the grouping patterns are study-specific.
2. **Linux — 01_06, 01_07, 01_08:** assign taxonomy and cluster OTUs.
   Supply the SILVA/UNITE directory as the first argument to 01_06 (or set
   `REFERENCE_DB_DIR`) and the Claident database **prefix** as the first
   argument to 01_07 (or set `CLAIDENT_DB`). Required reference releases are
   listed in the preprocessing environment record below.
   VSEARCH must be on `PATH`, or specified by `VSEARCH_PATH` for 01_08.
3. **Mac — 01_09:** normalize sample IDs and copy filtered tables from
   `Output/01_Data_processing/<dataset>/` to `Data/<dataset>/Seqdata/`.
   Also transfer the taxonomy and representative FASTA files written to
   `Data/` by 01_08. Continue at **Start from processed data**, step 1.

## Analysis environments

R, package, external-software and reference-database versions are recorded in:

- [Version_info_preprocessing.txt](Version_info/Version_info_preprocessing.txt)
  — Linux: 01_01–01_03 and 01_06–01_08.
- [Version_info_main_analysis.txt](Version_info/Version_info_main_analysis.txt)
  — Mac: 01_04–01_05, 01_09–01_10 and stages 02–08.


## Repository layout and execution

- `Scripts/01_*`–`08_*`: numbered workflow stages shown above.
- `Function/`: helper functions sourced by the analysis scripts.
- `Data/`: supplied inputs; `Output/`: generated results (not tracked on GitHub).
- `version_info/`: environment records.

Except for run-level 01_01–01_03, execute scripts from the repository root.
There is no single master script: datasets and analysis combinations must be
selected explicitly, and root annotation and tree reconstruction include
manual steps.

## Code provenance and citation

Some run-level sequence-processing scripts were developed by Hiroaki Fujita
and adapted by Akifumi Murata. Retain the original headers and attribution
when redistributing code.

Please cite: Murata et al. *Environmental filtering and host identity
collectively shape root-associated microbiomes of Ericaceae and
ectomycorrhizal plants in fumarole fields.*


This repository contains the scripts used for the study **Environmental
filtering and host identity collectively shape root-associated microbiomes of
Ericaceae and ectomycorrhizal plants in fumarole fields**.

The workflow covers sequence processing, plant-root identification, microbial
community filtering and characterization, soil analysis, community-level
tests, preference analysis, and fungal phylogenetic analysis.

## Data availability

Raw sequencing reads are not stored in this repository. They are available
from the DDBJ BioProject **PRJDB42913**. Processed sequence tables, OTU
taxonomy tables, and reference databases are also excluded from GitHub.

The sequence-processing scripts document how the processed tables were
generated from demultiplexed FASTQ files. Analyses that start from processed
sequence or taxonomy data require users to regenerate those files or obtain
access under the study's data-availability conditions.

The following external taxonomy databases are required but not distributed:

- SILVA `silva_nr99_v138.2_toSpecies_trainset.fa.gz`
- UNITE `sh_general_release_dynamic_s_19.02.2025.fasta`
- Claident `overall_genus` database

Database locations are supplied locally as described in
`Scripts/01_Data_processing/README.md`.

## Repository structure

```text
.
├── Scripts/
│   ├── 01_Data_processing/
│   ├── 02_Plant_root_identification/
│   ├── 03_Data_filtering/
│   ├── 04_Soil_analysis/
│   ├── 05_Community_characterization/
│   ├── 06_Community_analysis/
│   ├── 07_Preference_analysis/
│   └── 08_Phylogenetic_analysis/
├── Function/                  # helper functions sourced by analysis scripts
├── Data/                      # local inputs; not tracked on GitHub
├── Output/                    # generated results; not tracked on GitHub
└── Solfatara_Analysis.Rproj
```

Run scripts from the repository root so that paths created with the `here` R
package resolve correctly. The scripts create most output directories when
needed. Files in `Function/` are sourced by the analysis scripts and are not
run independently.

## Workflow

### 1 Sequence-data processing

Run the first three steps separately for each sequencing run and marker:

```text
01_01_cutadapt_R.sh
  -> 01_02_readQC_dada2.sh
  -> 01_03_Denoising_to_prechimera.sh
  -> stall_no_rmchimera.rds
```

The remaining steps combine runs before contaminant removal and downstream
processing:

```text
01_04_Merge_runs.R
  -> merge runs and remove chimeras
01_05_Decontamination.R
  -> prevalence-based decontamination for Prokaryote and Fungi
01_06_Taxonomy_annotation.R
  -> prepare ASV files and assign Prokaryote and Fungi taxonomy with DADA2
01_07_Claident_plant_annotation.sh
  -> assign Plant taxonomy with Claident
01_08_OTU_clustering.R
  -> cluster ASVs at 97% identity and retain representative taxonomy
01_10_Converting_prefix.R
  -> synchronize IDs in the final OTU table, taxonomy, and centroid FASTA
01_09_Coverage_rarefaction.R
  -> coverage-based rarefaction for Prokaryote and Fungi
```

All ASV and OTU identifiers remain `X_` through OTU clustering. Before
rarefaction and downstream analysis, `01_10_Converting_prefix.R` changes
Prokaryote IDs to `P_` and Fungi IDs to `F_`; Plant IDs remain `X_`.

The run-level pre-chimera tables and other large intermediate files are not
distributed. Their expected local paths and marker-specific parameters are
documented in `Scripts/01_Data_processing/README.md` and
`parameterList.tsv`.

### 2 Plant-root identification

Scripts in `02_Plant_root_identification/` identify plant OTUs in root samples
using reference leaf samples:

1. Prepare root and reference-leaf sequence tables and FASTA files.
2. Run local BLAST with the leaf sequences as the reference database.
3. Combine BLAST matches with leaf taxonomy and supplementary top-hit
   annotations.
4. Retain confidently assigned hosts and generate metadata for microbial
   analyses.
5. Plot host composition across sites and habitats.

The Plant input `seqOTUtab_filtered.rds` is filtered only to retain
Viridiplantae OTUs during sequence processing. Removal of negative controls
and selection among `_G`, `_G2`, `_G3`, and `_G4` gleaning replicates are
performed in `02_01_Process_seq_and_taxa.R`. This keeps taxonomy filtering
separate from sample-level processing.

`02_03_Root_taxonomy_annotation.R` reads the manually curated file
`LBlast_noref_ano.csv` for root OTUs without a reference-leaf match. This
manual annotation step must be completed before continuing to host filtering.

### 3 Community-data filtering

Scripts in `03_Data_filtering/` merge root and soil tables where required,
match sequence tables to sample metadata, and filter OTUs by occurrence. The
analyses use thresholds of 1, 3, or 5 samples. Set `data_type`, `sample_type`,
and `th` at the top of the relevant script and rerun it for each required
combination.

### 4 Soil analysis

Scripts in `04_Soil_analysis/` perform blank correction and limit-of-detection
filtering of exchangeable cations, analyze soil pH, compare habitats, and run
PCA, PERMANOVA, and multivariate-dispersion tests. The final scripts in this
directory format the corresponding supplementary tables as Word documents.

Run the processing scripts before their visualization and table scripts:

```text
04_01 and 04_02 -> 04_03 and 04_04 -> 04_05 -> 04_06 to 04_08
```

### 5 Community characterization

Scripts in `05_Community_characterization/` summarize microbial taxonomic
composition, estimate alpha diversity, test habitat and sample-type effects,
and create the corresponding figures. Set `data_type`, `sample_type`, and the
taxonomic rank where indicated before each run.

### 6 Community analysis

`06_00_Setup.R` loads a filtered community table and its metadata and computes
the Sorensen dissimilarity matrix. The other scripts use this setup to run
PERMANOVA, variation partitioning, PERMDISP, and PCoA, followed by
supplementary-table generation.

Set `data_type` and `sample_type` in `06_00_Setup.R`, then run the applicable
scripts in numerical order. Table scripts should be run only after their
corresponding analysis results have been generated.

### 7 Preference analysis

`07_00_Setup.R` defines the dataset, filtering threshold, preference method,
focal variable, and interaction direction. The remaining scripts calculate
two-dimensional preference and d-prime statistics, combine habitat and host
preference estimates, generate figures, test their correlation, and evaluate
robustness across OTU-filtering thresholds.

The main settings are:

- `data_type`: `Prokaryote` or `Fungi`
- `analysis`: `2DP` or `dprime`
- `focus`: `habitat` or `host` for 2DP
- `direction`: `microbe` or `host` for d-prime
- `threshold`: `1`, `3`, or `5`

Rerun the calculation scripts for every combination required by the figure or
robustness script.

### 8 Phylogenetic analysis

The fungal phylogenetic workflow is:

```text
08_01_Select_sequences.R
  -> 08_02_Sequence_alignment.sh       (MAFFT and trimAl)
  -> 08_03_MEGA12_NJ_tree_settings.txt (manual tree reconstruction in MEGA)
  -> 08_04_Phylotree_visualization.R
  -> 08_05_Phylogenetic_signal.R
```

Hyaloscyphaceae OTUs are selected and combined with reference sequences,
aligned with MAFFT, and trimmed with trimAl at a gap threshold of 0.6. The
Neighbor-Joining tree is reconstructed in MEGA 12.1 using the Kimura
two-parameter model, pairwise deletion, and 1,000 bootstrap replicates. Export
the tree as `Newick_Pairwise_0.6.nwk` before running the final two R scripts.

## Software requirements

The analysis scripts were developed with R 4.5.3. Major R dependencies include
`ape`, `Biostrings`, `bipartite`, `broom`, `dada2`, `decontam`, `DHARMa`,
`doParallel`, `doRNG`, `effectsize`, `emmeans`, `flextable`, `foreach`,
`ggh4x`, `ggnewscale`, `ggplot2`, `ggpubr`, `ggrepel`, `ggstar`, `ggtext`,
`ggtree`, `glmmTMB`, `here`, `officer`, `patchwork`, `phytools`, `rstatix`,
`seqinr`, `tidyverse`, `treeio`, and `vegan`.

External command-line programs include Cutadapt, FastQC, MultiQC, GNU
Parallel, VSEARCH, Claident, BLAST+, MAFFT 7.526, and trimAl 1.5.1. MEGA 12.1
is used interactively for phylogenetic-tree reconstruction.

For exact reproducibility, record package and system versions after creating
the analysis environment, for example with `sessionInfo()` in R and the
version commands provided by the external programs.

## Running the analysis

1. Clone the repository and open `Solfatara_Analysis.Rproj` or set the working
   directory to the repository root.
2. Install the required R packages and external programs.
3. Obtain the required input data and reference databases and place them at
   the paths documented by the scripts.
4. Edit only the settings block near the top of each script for the desired
   dataset and analysis combination.
5. Run scripts in the stage order described above. Generated files are written
   below `Output/`.

Example:

```bash
Rscript Scripts/06_Community_analysis/06_01_PERMANOVA.R
```

The repository does not provide a single master script because several stages
must be rerun for multiple marker, sample-type, and model combinations, and
the BLAST and MEGA stages include documented manual operations.

## Code provenance

Some run-level sequence-processing scripts were originally developed by
Hiroaki fujita and adapted for this study by Akifumi Murata. Replace the
placeholder with the preferred author attribution before release and retain
the original file headers.

## Citation and license

If you use this workflow, please cite:

> Murata et al. Environmental filtering and host identity collectively shape
> root-associated microbiomes of Ericaceae and ectomycorrhizal plants in
> fumarole fields. 

Add a repository license only after confirming that all code contributors
agree to distribute their contributions under that license.

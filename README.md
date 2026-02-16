# Project README

This repository contains analysis scripts and supporting files for genomic analyses, including centromere identification, Fst analysis, NLR characterisation, PCA workflows, pangenome graph construction, pangrowth analysis, and whole-genome alignments.

---

## Repository Structure

```
Code/
└── SuppScripts/
    ├── Centromeres/
    │   ├── MUMMER_alignments.txt
    │   ├── ProcessTRASHandMUMMER20250416.R
    │   ├── RunTRASHall.sh
    │   └── runCentIERall.sh
    │
    ├── Fst/
    │   ├── FstAnalysis.Rmd
    │   ├── MN106Ref_GenomicClasses.txt
    │   ├── MN106_GenePoorRegions_sorted.bed
    │   ├── MN106_SubTeloGenicRegions_sorted.bed
    │   ├── Tarvense_GBS_WGS_merged_GenePoorGDS
    │   ├── Tarvense_GBS_WGS_merged_GeneRichGDS
    │   └── well_and_origin.csv
    │
    ├── NLRs/
    │   ├── Full_NLR_analysis_code20250429.R
    │   ├── GOAnalysisNLROGs20250904.R
    │   ├── NLR_bash_code_20250429Final.txt
    │   ├── Pennycress_genespace.R
    │   └── combBed.tar.gz
    │
    ├── PCA/
    │   ├── Pennycress_PCA.R
    │   ├── Preprocess_genotypes.txt
    │   ├── README.txt
    │   ├── align_penncressGBS.sh
    │   ├── align_polishing_libs20241101.sh
    │   ├── call_polishing_libraries_20241105.sh
    │   ├── callsites_pennycressGBS.sh
    │   ├── demultiplex_pennycressGBS.sh
    │   └── filter_pennycressGBS.sh
    │
    ├── Pangenome_graph/
    │   ├── PG_01_minigraph-cactus_Chr01.sh
    │   └── PG_02_vg_combine_and_haplos.sh
    │
    ├── Pangrowth/
    │   ├── KMER_01_pangrowth.sh
    │   └── KMER_02_pangrowth_plot.R
    │
    └── WGAs/
        ├── WGA_01_minimap2_wga.sh
        ├── WGA_02_pafr_viz_and_write_filt.R
        ├── WGA_03_syri.sh
        └── WGA_04_syri_output_viz.R
```

---

## Directory Descriptions

### Centromeres

Scripts and outputs for identifying and analysing putative centromeres.

* MUMMER alignments for structural analysis
* TRASH and CentIER pipelines
* Downstream processing and plotting

### Fst

Population differentiation and genomic region classification analyses.

### NLRs

NLR gene identification, entropy analysis, and GO enrichment.

### PCA

Variant processing and principal component analysis workflows.

### Pangenome Graph

Construction and processing of chromosome-level pangenome graphs.

### Pangrowth

K-mer based pangrowth analysis and visualisation.

### WGAs

Whole-genome alignment pipelines and visualisation.

---

## Notes

* Scripts are grouped by analysis type.
* Some directories contain intermediate or processed data files required for downstream steps.
* See subdirectory README files where available for pipeline details.

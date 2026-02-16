README

Code/SuppScripts
    Centromeres
        MUMMER_alignments.txt - alignments of centromeric regions to analyze structure
        ProcessTRASHandMUMMER20250416.R - downstream analysis and plotting of TRASH and CentIER output
        RunTRASHall.sh - code to identify putative centromeres with TRASH
        runCentIERall.sh - code to identify putative centromeres with CentIER
Fst
    FstAnalysis.Rmd
    MN106Ref_GenomicClasses.txt
    MN106_GenePoorRegions_sorted.bed
    MN106_SubTeloGenicRegions_sorted.bed
    Tarvense_GBS_WGS_merged_GenePoorGDS
    Tarvense_GBS_WGS_merged_GeneRichGDS
    well_and_origin.csv
NLRs
    Full_NLR_analysis_code20250429.R - code to identify and process NLRs, including entropy analysis
    GOAnalysisNLROGs20250904.R - GO analysis of terms enriched in the dispensible genome
    NLR_bash_code_20250429Final.txt - code for hmmer identification of NLRs
    Pennycress_genespace.R - code to run GENESPACE on Thlaspi and other accessions
    combBed.tar.gz - GENESPACE combbed file used in NLR analyses (including T. arvense accessions and TAIR11)
PCA
    Pennycress_PCA.R - code to generate PCA with T. arvense accessions
    Preprocess_genotypes.txt - code to prepare genotypes for PCA
    README.txt - alignment and variant calling pipeline 
    align_penncressGBS.sh
    align_polishing_libs20241101.sh
    call_polishing_libraries_20241105.sh
    callsites_pennycressGBS.sh
    demultiplex_pennycressGBS.sh
    filter_pennycressGBS.sh
Pangenome_graph
    PG_01_minigraph-cactus_Chr01.sh
    PG_02_vg_combine_and_haplos.sh
Pangrowth
    KMER_01_pangrowth.sh
    KMER_02_pangrowth_plot.R
WGAs
    WGA_01_minimap2_wga.sh
    WGA_02_pafr_viz_and_write_filt.R
    WGA_03_syri.sh
    WGA_04_syri_output_viz.R
###  Script by P. Grabowski, edited by A. Harder  ###
###     Visualization of minimap2 WGA results     ###

setwd('~/mm2_output/')

### LOAD PACKAGES ###
# install.packages('pafr')
# install.packages('ggpubr')
library(pafr)
library(ggplot2)
library(ggpubr)
library(scales)

## set filtering criteria to test
id.cut <- 0.80
len.cut <- 1e3
min.mapq <- 2

### SET INPUTS ###
pafs <- list.files(path = '.', pattern = "*-f_0.paf")
for(p in pafs){
  print(p)
  paf_file <- p
  samp1 <- unlist(strsplit(paf_file, split = '_'))[1]
  samp2 <- unlist(strsplit(paf_file, split = '_'))[3]
  mm2 <- gsub('.paf', '', unlist(strsplit(paf_file, split = '_'))[4])
  ali_0 <- read_paf(paf_file)
  
  ### SET VARIABLES ###
  ref_chrs <- sort(unique(ali_0$tname))
  quer_chrs <- sort(unique(ali_0$qname))
  
  # filtering parameters 
  per_idy_cut <- id.cut
  len_cut <- len.cut
  
  ### SET OUTPUTS ###
  out_dir <- '../R_figures/' 
  filt_dir <- '../filtered_mm2_output/'
  plot_pre <- gsub('.paf', '', paf_file)
  
  # adjust names as need be
  filt.paf <- paste0(filt_dir,'filtered_',plot_pre,'_',len.cut,'_',id.cut,'_',min.mapq,'.paf')
  unfilt_dot_pdf <- paste0(out_dir, '01_', plot_pre, '_unfilt_dotplot.pdf') 
  filt_dot_pdf <- paste0(out_dir, '02_', plot_pre, '_filt_dotplot.pdf')
  unfilt_cov_pdf <- paste0(out_dir, '03_', plot_pre, '_unfilt_align_coverage.pdf') 
  filt_cov_pdf <- paste0(out_dir, '04_', plot_pre, '_filt_align_coverage.pdf')
  perc.id.pdf <- paste0(out_dir, '05_', plot_pre, '_percent_id_freq.pdf')
  
  # add percent identity info
  ali_0$PER_IDY <- ali_0$nmatch / ali_0$alen
  
  # optional: remove scaffolds
  # replace ali_chr_0 with ali_0 below if don't do this step 
  ali_chr_0 <- ali_0[which(ali_0$qname %in% quer_chrs & ali_0$tname %in% ref_chrs), ]
  
  # unfiltered dotplot
  unfilt_keep <- list(intersect(quer_chrs, unique(ali_chr_0$qname)), intersect(ref_chrs, unique(ali_chr_0$tname)))
  unfilt_dot <- dotplot(ali_chr_0, label_seqs = T, order_by = 'provided', ordering = unfilt_keep, xlab = samp2, ylab = samp1)
  pdf(unfilt_dot_pdf, width = 25, height = 25)
  print(unfilt_dot)
  dev.off()
  
  # unfiltered coverage plot
  unfilt_cov <- plot_coverage(ali_chr_0, target = T)
  pdf(unfilt_cov_pdf)
  print(unfilt_cov)
  dev.off()
  
  # Look at filtered results
  # optional: remove secondary alignments
  prim_ali <- filter_secondary_alignments(ali_chr_0)
  ali_filt <- prim_ali[which(prim_ali$PER_IDY > per_idy_cut & prim_ali$alen > len_cut & prim_ali$mapq >= 2), ]
  
  hist(ali_filt$mapq, xlim = c(0, 57), breaks = 100, ylim = c(0, 100), main = samp2)
  ## read in raw PAF file to keep all formatting, subset based on row #s
  raw.paf <- read.table(paf_file, sep = '?')
  to.write <- raw.paf[rownames(ali_filt),]
  write.table(to.write, filt.paf, col.names = FALSE, row.names = FALSE, quote = FALSE)
  
  # dot plot
  to_keep <- list(intersect(quer_chrs, unique(ali_filt$qname)),
                  intersect(ref_chrs, unique(ali_filt$tname)))
  tot_dot <- dotplot(ali_filt, label_seqs = T, order_by = 'provided', ordering = to_keep, xlab = samp2, ylab = samp1)
  pdf(filt_dot_pdf, width = 25, height = 25) 
  print(tot_dot)
  dev.off()
  
  # coverage plot
  good_cov <- plot_coverage(ali_filt, target = T) 
  pdf(filt_cov_pdf)
  print(good_cov)
  dev.off()
}
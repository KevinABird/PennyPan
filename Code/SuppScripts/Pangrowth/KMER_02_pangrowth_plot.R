setwd('~/pangrowth_output/')

fns <- list.files(path = 'NLR_specific/', pattern = '_nonhvNLRs_growth.txt')
ks <- gsub('k', '', do.call(rbind, strsplit(fns, split = '_', fixed = TRUE))[,3])

##### Plot growth (fitted values) data #####
method <- 'growth'

for(k in ks){
  wg <- t(read.table(paste0('wg_genes_NLRs/pennycress_n7_k',k,'_growth.txt')))
  idx <- seq(1:length(wg))
  genes <- t(read.table(paste0('wg_genes_NLRs/pennycress_n7_k',k,'_genes_growth.txt')))
  NLRs <- t(read.table(paste0('wg_genes_NLRs/pennycress_n7_k',k,'_NLR_growth.txt')))
  hvNLRs <- t(read.table(paste0('NLR_specific/pennycress_n7_k',k,'_hvNLRs_growth.txt')))
  
  wg.norm <- wg/(max(wg))
  genes.norm <- genes/max(genes)
  NLRs.norm <- NLRs/max(NLRs)
  hvNLRs.norm <- hvNLRs/(max(hvNLRs))
  
  wg.col <- 'dodgerblue2'
  genes.col <- 'darkorange2'
  NLRs.col <- 'darkgoldenrod3'
  hvNLRs.col <- 'goldenrod1'
  
  pdf(paste0('../R_figures/pennycress_kmer-combined_growth-k',k,'.pdf'), width = 6, height = 6)
  
  plot(idx, wg.norm, ylim = c(0,1), xlab = 'Number of haplotypes in collection',
       ylab = 'Proportional k-mer distribution',
       col = 'transparent')
  
  points(idx, wg.norm, pch = 16, col = wg.col)
  lines(idx, wg.norm, col = wg.col)
  
  points(idx, genes.norm, pch = 16, col = genes.col)
  lines(idx, genes.norm, col = genes.col)
  
  points(idx, NLRs.norm, pch = 16, col = NLRs.col)
  lines(idx, NLRs.norm, col = NLRs.col)
  
  points(idx, hvNLRs.norm, pch = 16, col = hvNLRs.col)
  lines(idx, hvNLRs.norm, col = hvNLRs.col)
  
  legend('bottomright', legend = c('Whole genome', 'Genes', 'NLRs', 'hvNLRs'),
         pch = 16, lty = 1, col = c(wg.col, genes.col, NLRs.col, hvNLRs.col), bty = 'n',
         inset = 0.02)
  dev.off()  
}

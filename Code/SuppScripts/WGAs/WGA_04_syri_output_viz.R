##### Quantify interchromosomal alignments and unaligned segments for full genome alignments #####
setwd('~/syri_output/')
library(ggplot2)
library(paletteer)

## total reference genome length
chrom.lens <- read.table('../background/Tarvensevar_MN106_872_v4.0.fa.fai')
colnames(chrom.lens) <- c('chrom','length')
tot.ref.len <- sum(chrom.lens$length)

fns <- list.files()

plot.chrom <- 'all chroms'
## turn off this loop if you wanna assess all chroms simultaneously
for(plot.chrom in chrom.lens$chrom){

  OUT <- NULL
  OUT1 <- NULL ## set up stacked bar plot data
  for(f in fns){
    print(f)
    dat <- read.table(f)
    samp2 <- unlist(strsplit(f, split = '_'))[1]
    ## positions all 1-based
    colnames(dat) <- c('ref.chrom','ref.s','ref.e','ref.allele',
                       'comp.allele','comp.chrom','comp.s','comp.e',
                       'uniq.id','parent.id','annot','copy.status')
    for(c in c(2,3,7,8)){
      dat[,c] <- suppressWarnings(as.numeric(dat[,c]))
    } ## NAs introduced, but they're not NOTAL segments
    
    dat$ref.len <- abs(dat$ref.e - dat$ref.s) + 1
    dat$comp.len <- abs(dat$comp.e - dat$comp.s) + 1
    dat[is.na(dat$ref.len), 'ref.len'] <- 0
    dat[is.na(dat$comp.len), 'comp.len'] <- 0
    
    ## just plot one chromosome or turn off
    # dat <- dat[dat$comp.chrom == plot.chrom & dat$ref.chrom %in% c(plot.chrom, '-'),]
    
    dat <- dat[dat$comp.chrom %in% c(plot.chrom, '-') & dat$ref.chrom %in% c(plot.chrom, '-'),]
    
    ## save data for stacked bar plot
    for(cat in c('SYNAL','HDR','DUPAL','INVAL','INVDPAL','INVTRAL','NOTAL','INS','CPG')){
      save <- c(samp2, cat, sum(dat[dat$annot == cat, 'comp.len']))
      OUT1 <- rbind(OUT1, save)
    }
    
    ## keep regions & alignments (subsets of regions), not SNPs, etc.
    regions <- dat[dat$annot %in% c('SYN','INV','TRANS','INVTR','DUP','INVDP'),]
    syn.sum <- sum(dat[dat$annot == 'SYN', 'ref.len'])
    hdr.sum <- sum(dat[dat$annot == 'HDR', 'ref.len'])
    synal.sum <- sum(dat[dat$annot == 'SYNAL', 'ref.len'])
    notal.sum <- sum(dat[dat$annot == 'NOTAL', 'ref.len'], na.rm = TRUE)
    # print(paste0(samp2,' - SYN prop: ',round(syn.sum/tot.ref.len, digits = 2),
    #              ' ; HDR/SYN: ',round(hdr.sum/syn.sum, digits = 2),
    #              ' ; SYNAL/SYN: ',round(synal.sum/syn.sum, digits = 2),
    #              ' ; NOTAL prop: ',round(notal.sum/tot.ref.len, digits = 2)))
    print(paste0(samp2,' - ',hdr.sum/tot.ref.len))
    ref.region.cov <- sum(regions$ref.len)/tot.ref.len
    aligns <- dat[dat$annot %in% c('SYNAL','INVAL','TRANSAL','INVTRAL','DUPAL','INVDPAL'),]
    ref.align.cov <- sum(aligns$ref.len)/tot.ref.len
    
    ## split out unaligned segments and interchromosomal alignments
    notal <- dat[dat$annot == 'NOTAL',]
    ic <- aligns[aligns$ref.chrom != aligns$comp.chrom,]
  
    ## calculate total unaligned length
    full.notal.ref <- sum(notal$ref.len, na.rm = TRUE)/tot.ref.len
    full.notal.comp <- sum((notal$comp.e - notal$comp.s + 1), na.rm = TRUE)/tot.ref.len
    notal.prop <- sum(notal$ref.len, na.rm = TRUE)/tot.ref.len
    
    ## calculate total interchromosomal alignments
    ic.prop <- sum(ic$ref.len, na.rm = TRUE)/tot.ref.len
    
    save <- c(samp2, 'full.genome.align', ref.region.cov, ref.align.cov, 
              notal.prop, ic.prop)
    OUT <- rbind(OUT, save)
  }
  full <- as.data.frame(OUT)
  colnames(full) <- c('samp2', 'minimap.level', 'ref.region.cov', 'ref.align.cov', 
                      'notal.prop','ic.prop')
  for(c in 3:ncol(full)){
    full[,c] <- as.numeric(full[,c])
  }
  samps <- unique(full$samp2)
  
  bardat <- as.data.frame(OUT1)
  colnames(bardat) <- c('samp','alignment.category','length')
  bardat$length <- as.numeric(bardat$length)
  
  ## plot stacked bar plot
  ## choose categories to plot
  tmp <- bardat[bardat$alignment.category %in% c('SYNAL','HDR','DUPAL','INVAL','INVTRAL','NOTAL'),]
  
  gg <- ggplot(tmp, aes(fill = alignment.category, x = length, y = samp)) +
    geom_bar(position = 'stack', stat = 'identity') + 
    scale_x_continuous(labels=function(x)x/1e6) +
    scale_fill_paletteer_d("MetBrewer::Archambault") +
    ggtitle(paste0('SyRI results - ',plot.chrom)) +
    xlab('Length (Mb)') +
    ylab('')

  pdf(paste0('../R_figures/syri_alignment_category_length_totals_',plot.chrom,'.pdf'), width = 9, height = 5)
    gg
  dev.off()
}
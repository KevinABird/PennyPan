library(SNPRelate)
library(SeqArray)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)

setwd("~/Pennycress/polishing_library_SNP_calls")
vcf.fn <- "small_set_GBS_and_WGS_pennycress_SNPs.vcf"
# Reformat
snpgdsVCF2GDS(vcf.fn, "pennycress.gds", method="biallelic.only")
genofile <- snpgdsOpen("pennycress.gds")
select_samples<-read.gdsn(index.gdsn(genofile, "sample.id"))
select_samples<-select_samples[select_samples!="H08"] #remove blank sample
chr_list<-read.gdsn(index.gdsn(genofile, "snp.chromosome"))
chr_06_positions<-which(chr_list == "06")
read.gdsn(index.gdsn(genofile, "snp.id"))
set.seed(1000)



snpset <- snpgdsLDpruning(genofile, ld.threshold=.95, maf = 0.1, missing.rate = 0.5, sample.id = select_samples)

snpset_chr6 <- snpgdsLDpruning(genofile, ld.threshold=.95, maf = 0.1, missing.rate = 0.5, sample.id = select_samples,  snp.id = chr_06_positions)

snpset.id <- unlist(unname(snpset))
snpset.id.chr6 <- unlist(unname(snpset_chr6))

pca <- snpgdsPCA(genofile, snp.id=snpset.id, num.thread=2, sample.id = select_samples)
pc.percent <- pca$varprop*100
head(round(pc.percent, 2))
tab <- data.frame(sample.id = pca$sample.id,
                  EV1 = pca$eigenvect[,1],    # the first eigenvector
                  EV2 = pca$eigenvect[,2],    # the second eigenvector
                  EV3 = pca$eigenvect[,3],
                  stringsAsFactors = FALSE)

tab$is.ref<-rep(c(0,1), c(7, nrow(tab)-7))
tab$plot_label<-fcase(tab$is.ref==0, tab$sample.id)
tab$Armenia<-(tab$sample.id %in% c("E08", "F08", "G08", "H12", "A09", "B09", "C09"))
tab$is.ref.Armenia<-(fcase(tab$Armenia==T, 2,
                           tab$Armenia==F, tab$is.ref))



pca.chr6 <- snpgdsPCA(genofile, snp.id=snpset.id.chr6, num.thread=2, sample.id = select_samples)
pc.percent.chr6 <- pca.chr6$varprop*100
head(round(pc.percent.chr6, 2))
tab.chr6 <- data.frame(sample.id = pca.chr6$sample.id,
                  EV1 = pca.chr6$eigenvect[,1],    # the first eigenvector
                  EV2 = pca.chr6$eigenvect[,2],    # the second eigenvector
                  EV3 = pca.chr6$eigenvect[,3],
                  stringsAsFactors = FALSE)

tab.chr6$is.ref<-rep(c(0,1), c(7, nrow(tab.chr6)-7))
tab.chr6$plot_label<-fcase(tab.chr6$is.ref==0, tab.chr6$sample.id)
tab.chr6$Armenia<-(tab.chr6$sample.id %in% c("E08", "F08", "G08", "H12", "A09", "B09", "C09"))
tab.chr6$is.ref.Armenia<-(fcase(tab.chr6$Armenia==T, 2,
                                tab.chr6$Armenia==F, tab.chr6$is.ref))


pdf("SNPRelate_PCAs_20241126.pdf")
ggplot(tab, aes(x=EV1, y=EV2, colour = as.factor(is.ref), shape=as.factor(is.ref))) + 
  geom_hline(yintercept = 0, linewidth = .3)+
  geom_vline(xintercept = 0, linewidth = .3)+
  geom_point()+
  theme_light()+
  xlab("PC 1 (18.94%)")+
  ylab("PC 2 (8.45%)") +
  theme(legend.position="none") +
  geom_label_repel(label=tab$plot_label) +
  ggtitle("A. PCA of Thlaspi arvense diversity and sequenced reference genomes")+
  scale_color_brewer( palette = "Dark2") +
  scale_shape_manual(values=c(8,19))

ggplot(tab, aes(x=EV2, y=EV3, colour = as.factor(is.ref), shape=as.factor(is.ref))) + 
  geom_hline(yintercept = 0, linewidth = .3)+
  geom_vline(xintercept = 0, linewidth = .3)+
  geom_point()+
  theme_light()+
  xlab("PC 2 (8.45%)")+
  ylab("PC 3 (7.30%)") +
  theme(legend.position="none") +
  geom_label_repel(label=tab$plot_label) +
  ggtitle("B. PCA of Thlaspi arvense diversity and sequenced reference genomes")+
  scale_color_brewer( palette = "Dark2") +
  scale_shape_manual(values=c(8,19))


ggplot(tab, aes(x=EV1, y=EV2, colour = as.factor(is.ref.Armenia), shape=as.factor(is.ref))) + 
  geom_hline(yintercept = 0, linewidth = .3)+
  geom_vline(xintercept = 0, linewidth = .3)+
  geom_point()+
  theme_light()+
  xlab("PC 1 (18.94%)")+
  ylab("PC 2 (8.45%)") +
  theme(legend.position="none") +
  geom_label_repel(label=tab$plot_label) +
  ggtitle("C. PCA of Thlaspi arvense diversity and sequenced reference genomes highlighting Armenian samples")+
  scale_color_brewer( palette = "Dark2") +
  scale_shape_manual(values=c(8,19))


ggplot(tab.chr6, aes(x=EV1, y=EV2, colour = as.factor(is.ref.Armenia), shape=as.factor(is.ref))) + 
  geom_hline(yintercept = 0, linewidth = .3)+
  geom_vline(xintercept = 0, linewidth = .3)+
  geom_point()+
  theme_light()+
  xlab("PC 1 (18.94%)")+
  ylab("PC 2 (8.45%)") +
  theme(legend.position="none") +
  geom_label_repel(label=tab.chr6$plot_label) +
  ggtitle("D. PCA of Thlaspi arvense diversity and sequenced reference genomes on chromosome 6 only")+
  scale_color_brewer( palette = "Dark2") +
  scale_shape_manual(values=c(8,19))

dev.off()

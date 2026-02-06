##### Set working directory and load packages #####
setwd("~/NLRs")
library(data.table)
library(stringr)
library(rtracklayer)
library(GENESPACE)
library(GenomicRanges)
library(BiocGenerics)
library(Biostrings)
library(entropy)
library(tidyverse)
library(pwalign)
mafft="/usr/local/bin/mafft"
sessionInfo()



#R version 4.4.2 (2024-10-13)
#GENESPACE_1.3.1
#updated 20250213 to add A. lyrata and B. rapa
#sessionInfo()
#  [1] lubridate_1.9.4      forcats_1.0.0        dplyr_1.1.4          purrr_1.0.4          readr_2.1.5          tidyr_1.3.1          tibble_3.2.1        
#[8] ggplot2_3.5.1        tidyverse_2.0.0      entropy_1.3.1        Biostrings_2.72.1    XVector_0.44.0       GENESPACE_1.3.1      rtracklayer_1.64.0  
#[15] GenomicRanges_1.56.2 GenomeInfoDb_1.40.1  IRanges_2.38.1       S4Vectors_0.42.1     BiocGenerics_0.50.0  stringr_1.5.1        data.table_1.17.0   



####### IDENTIFY NLRS #######

##### Filter HMMER candidates #####

#see NLR_bash_code.txt to generate input file

#Collect hmmer hits, filter by e value, and output lists of sequence names
for (seqname in c("AK34W", "MN106", "MN134", "PI650286", "Ames", "Loretto", "Tibet")) {
  print(seqname)
  hmmer_out_name<-paste("~/NLRs/", seqname, "_NLR.out", sep="" ) #update path to match local path
  print(hmmer_out_name)
  hmmer_file<-fread(hmmer_out_name, skip = 14, fill=T, header = T)
  print(hmmer_file)
  hmmer_file$`E-value`<-as.numeric(hmmer_file$`E-value`)
  hmmer_file<-hmmer_file[`E-value`<0.01]
  outname<-paste(seqname, "_seqnames.txt", sep="")
  fwrite(as.list(hmmer_file$Sequence), outname, sep="\n")
}

for (seqname in c("TAIR", "Alyrata", "Brapa")) {
  print(seqname)
  hmmer_out_name<-paste(seqname, "_NLR.out", sep="" )
  print(hmmer_out_name)
  hmmer_file<-fread(hmmer_out_name, skip = 14, fill=T, header = T)
  print(hmmer_file)
  hmmer_file$`E-value`<-as.numeric(hmmer_file$`E-value`)
  hmmer_file<-hmmer_file[`E-value`<0.01]
  outname<-paste(seqname, "_seqnames.txt", sep="")
  fwrite(as.list(hmmer_file$Sequence), outname, sep="\n")
}
#This throws a warning but actually it's just stopping at the inclusion threshold


#see NLR_bash_code.txt to generate input file

for (seqname in c("AK34W", "MN106", "MN134", "PI650286", "Ames", "Loretto", "Tibet")) {
  print(seqname)
  BLASTname<-paste("BLAST_",seqname, "_full_hitdata.txt", sep="")
  BLAST<-fread(BLASTname, header = T, skip=7)
  BLAST<-BLAST[grepl(x=`Short name`, pattern="NB-ARC")]
  BLAST<-BLAST[Incomplete=="-"]
  BLAST[, gene_name:= tstrsplit(Query, " ", fixed = TRUE, keep = 3)]
  BLAST$gene_name<-gsub(">", "", BLAST$gene_name)
  BLAST[,transcript:=tstrsplit(gene_name, ".", fixed=TRUE, keep=3)]
  BLAST<-BLAST[transcript==1]
  BLAST[, pacid:= tstrsplit(Query, " ", fixed = TRUE, keep = 4)]
  BLAST$pacid<-gsub("pacid=", "", BLAST$pacid)
  GFFpath<-list.files(path="../genomes", pattern = seqname, full.names = T)
  GFFname<-list.files(path=GFFpath, pattern = ".gff3.gz", full.names = T)
  GFF<-as.data.table(readGFF(GFFname))
  GFF[type=="gene"]
  GFF<-GFF[pacid %in% BLAST$pacid & type == "mRNA"]
  print(table(GFF$seqid))
  outname<-paste(seqname, "_full_NLR_hits.gff", sep="")
  fwrite(GFF[,1:9],outname, sep="\t")
}


#slightly different code for Arabidopsis thaliana
seqname="TAIR"
BLASTname<-paste("BLAST_",seqname, "_full_hitdata.txt", sep="")
BLAST<-fread(BLASTname, header = T, skip=7)
BLAST<-BLAST[grepl(x=`Short name`, pattern="NB-ARC")]
BLAST<-BLAST[Incomplete=="-"]
BLAST[, gene_name:= tstrsplit(Query, " ", fixed = TRUE, keep = 3)]
BLAST$gene_name<-gsub(">", "", BLAST$gene_name)
BLAST[,transcript:=tstrsplit(gene_name, ".", fixed=TRUE, keep=2)]
BLAST<-BLAST[transcript==1]
BLAST[, pacid:= tstrsplit(Query, " ", fixed = TRUE, keep = 4)]
BLAST$pacid<-gsub("pacid=", "", BLAST$pacid)
GFFpath<-list.files(path="../genomes", pattern = seqname, full.names = T)
GFFname<-list.files(path=GFFpath, pattern = ".gff3.gz", full.names = T)
GFF<-as.data.table(readGFF(GFFname))
GFF[type=="gene"]
GFF<-GFF[pacid %in% BLAST$pacid & type == "mRNA"]
print(table(GFF$seqid))
outname<-paste(seqname, "_full_NLR_hits.gff", sep="")
fwrite(GFF[,1:9],outname, sep="\t")

#slightly different code for Arabidopsis lyrata
seqname="Alyrata"
BLASTname<-paste("BLAST_",seqname, "_full_hitdata.txt", sep="")
BLAST<-fread(BLASTname, header = T, skip=7)
BLAST<-BLAST[grepl(x=`Short name`, pattern="NB-ARC")]
BLAST<-BLAST[Incomplete=="-"]
BLAST[, gene_name:= tstrsplit(Query, " ", fixed = TRUE, keep = 3)]
BLAST$gene_name<-gsub(">", "", BLAST$gene_name)
BLAST[,transcript:=tstrsplit(gene_name, ".", fixed=TRUE, keep=2)]
BLAST<-BLAST[transcript=="t1"]
BLAST[, pacid:= tstrsplit(Query, " ", fixed = TRUE, keep = 4)]
BLAST$pacid<-gsub("pacid=", "", BLAST$pacid)
GFFpath<-list.files(path="../genomes", pattern = seqname, full.names = T)
GFFname<-list.files(path=GFFpath, pattern = ".gff3.gz", full.names = T)
GFF<-as.data.table(readGFF(GFFname))
GFF[type=="gene"]
GFF<-GFF[pacid %in% BLAST$pacid & type == "mRNA"]
print(table(GFF$seqid))
outname<-paste(seqname, "_full_NLR_hits.gff", sep="")
fwrite(GFF[,1:9],outname, sep="\t")


#slightly different code for Brassica rapa
seqname="Brapa"
BLASTname<-paste("BLAST_",seqname, "_full_hitdata.txt", sep="")
BLAST<-fread(BLASTname, header = T, skip=7)
BLAST<-BLAST[grepl(x=`Short name`, pattern="NB-ARC")]
BLAST<-BLAST[Incomplete=="-"]
BLAST[, gene_name:= tstrsplit(Query, " ", fixed = TRUE, keep = 3)]
BLAST$gene_name<-gsub(">", "", BLAST$gene_name)
BLAST[,transcript:=tstrsplit(gene_name, ".", fixed=TRUE, keep=3)]
BLAST<-BLAST[transcript=="1"]
BLAST[, pacid:= tstrsplit(Query, " ", fixed = TRUE, keep = 4)]
BLAST$pacid<-gsub("pacid=", "", BLAST$pacid)
GFFpath<-list.files(path="../genomes", pattern = seqname, full.names = T)
GFFname<-list.files(path=GFFpath, pattern = ".gff3.gz", full.names = T)
GFF<-as.data.table(readGFF(GFFname))
GFF[type=="gene"]
GFF<-GFF[pacid %in% BLAST$pacid & type == "mRNA"]
print(table(GFF$seqid))
outname<-paste(seqname, "_full_NLR_hits.gff", sep="")
fwrite(GFF[,1:9],outname, sep="\t")

##### NLR ORTHOGROUPS #####



##### Extract and count orthogroups containing pennycress NLRs #####

#import orthogroups from genespace
orthogroups<-GENESPACE::parse_hogs("DataOutputs/N0.tsv")
orthogroups[,Orthogroup:=tstrsplit(hogID, split=" ",fixed=TRUE, keep=2)]
setnames(orthogroups, "id", "ID")
#designate Arabidopsis hvNLRs
Arabidopsis_hvNLRs<-c("AT5G48620","AT5G46520","AT5G46510","AT5G43740","AT5G43470","AT5G41750","AT5G41740","AT5G38350","AT4G16950","AT4G16920","AT4G16890","AT4G16860","AT3G46530","AT3G44670","AT3G44630","AT3G44480","AT3G44400","AT2G14080","AT1G69550","AT1G62630","AT1G61310","AT1G61300","AT1G61190","AT1G61180","AT1G59218","AT1G59124","AT1G58848","AT1G58807","AT1G58602","AT1G31540") #identified from Prigozhin et al. 2021 https://doi.org/10.1093/plcell/koab013

#Collect all NLR hits from all taxa
NLR_gffs<-list.files(path="DataOutputs/", pattern="full_NLR_hits.gff", full.names = T)
NLR_gffs<-NLR_gffs[!str_detect(NLR_gffs,pattern="Brapa")]
NLR_gffs<-NLR_gffs[!str_detect(NLR_gffs,pattern="Alyrata")]
names(NLR_gffs)<-sapply(strsplit(NLR_gffs, split="\\_"),"[", 1)
names(NLR_gffs)<-sapply(strsplit(names(NLR_gffs), split="\\//"),"[", 2)
BigGff<-rbindlist(lapply(names(NLR_gffs), function(i)
  data.table(file_name=i,
             fread(file=NLR_gffs[i]))
))

#tidy gene ID names
BigGff$ID<-gsub(".1.p.v1.1", "", BigGff$ID)
BigGff$ID<-gsub(".1.v1.1", "", BigGff$ID)
BigGff$ID<-gsub(".1.v4.1", "", BigGff$ID)
BigGff$ID<-gsub(".1.Araport11.447", "", BigGff$ID)
BigGff$ID<-gsub(".t1.v2.1", "", BigGff$ID)
BigGff$ID<-gsub(".1.v2.1", "", BigGff$ID)

#select subset of orthogroups that HMMER and BLAST identified as containing NLRs 
HMMER_orthogroups<-orthogroups[ID %in% BigGff$ID]

#count NLRs in each genome
HMMER_orthogroups_counts<-as.data.table(as.data.frame.matrix(table(HMMER_orthogroups$Orthogroup, HMMER_orthogroups$genome)))
rownames(as.data.frame.matrix(table(HMMER_orthogroups$Orthogroup, HMMER_orthogroups$genome)))
HMMER_orthogroups_counts$Orthogroup<-rownames(as.data.frame.matrix(table(HMMER_orthogroups$Orthogroup, HMMER_orthogroups$genome)))

HMMER_orthogroups_counts[, variable_size_in_pennycress := apply(.SD, 1, function(row) !all(row == row[1])), .SDcols=(colnames(HMMER_orthogroups_counts))[2:8]]
HMMER_orthogroups_counts[, absent_in_pennycress := apply(.SD, 1, function(row) all(row == 0)), .SDcols=(colnames(HMMER_orthogroups_counts))[2:8]]

HMMER_orthogroups[, arabidopsis_hvNLR := c(0,1)[grepl(pattern = 
                                                        paste(Arabidopsis_hvNLRs, collapse="|"), ID) + 1L]]
HMMER_orthogroups_counts[, includes_arabidopsis_hvNLR := Orthogroup %in% HMMER_orthogroups[arabidopsis_hvNLR==1]$Orthogroup]
HMMER_orthogroups_counts[,size_in_pennycress:=rowSums(HMMER_orthogroups_counts[,2:8])]
fwrite(HMMER_orthogroups_counts, "DataOutputs/HMMER_orthogroup_counts_20250430.csv")

##### Add orthogroup IDs to gffs ##### 


GFF_with_orthogroups<-orthogroups[BigGff, on = "ID"]
rm(BigGff)
rm(orthogroups)

GFF_with_orthogroups[, includes_arabidopsis_hvNLR := Orthogroup %in% HMMER_orthogroups[arabidopsis_hvNLR==1]$Orthogroup]
#fwrite(GFF_with_orthogroups, "DataOutputs/TAIR_and_thlaspi_genes_by_orthogroup_and_hvNLR_20250430.csv")

##### Add NLR N-terminal domain class to GFF data #####
#GFF_with_orthogroups<-fread("DataOutputs/TAIR_and_thlaspi_genes_by_orthogroup_and_hvNLR_20250430.csv")
full_BLASTs<-list.files(path="DataOutputs", pattern="_full_hitdata.txt", full.names = T)
names(full_BLASTs)<-sapply(strsplit(full_BLASTs, split="\\_"),"[", 2)
all_full_BLAST<-rbindlist(fill=T, lapply(names(full_BLASTs), function(i)
  data.table(file_name=i,
             fread(file=full_BLASTs[i],, header = T, skip=7))
))
all_full_BLAST[, ID:= tstrsplit(Query, "ID=", fixed = TRUE, keep = 2)]
all_full_BLAST$ID<-tstrsplit(all_full_BLAST$ID, " ", fixed=TRUE, keep=1)
all_full_BLAST$ID<-gsub(".\\d.p.v1.1", "", all_full_BLAST$ID)
all_full_BLAST$ID<-gsub(".\\d.v4.1", "", all_full_BLAST$ID)
all_full_BLAST$ID<-gsub(".\\d.p.v4.1", "", all_full_BLAST$ID)
all_full_BLAST$ID<-gsub(".\\d.Araport11.447", "", all_full_BLAST$ID)
all_full_BLAST$ID<-gsub(".\\d.t1.v2.1", "", all_full_BLAST$ID)
all_full_BLAST$ID<-gsub(".\\d.\\d.p.v2.1", "", all_full_BLAST$ID)
all_full_BLAST[, TIR_domain:=fcase(
  grepl(x = `Short name`, pattern = "TIR"), 1,
  default=0)]
all_full_BLAST[, CC_domain:=fcase(
  grepl(x = `Definition`, pattern = "Coiled-coil"), 1,
  default=0)]
all_full_BLAST[, RPP8_domain:=fcase(
  grepl(x = `Definition`, pattern = "RPW8"), 1,
  default=0)]
TIR_NLRs<-unique(all_full_BLAST[TIR_domain==1,c("ID", "TIR_domain")])
CC_NLRs<-unique(all_full_BLAST[CC_domain==1,c("ID", "CC_domain")])
RPP8_NLRs<-unique(all_full_BLAST[RPP8_domain==1,c("ID", "RPP8_domain")])

GFF_with_orthogroups<-TIR_NLRs[GFF_with_orthogroups, on="ID"]
GFF_with_orthogroups<-CC_NLRs[GFF_with_orthogroups, on="ID"]
GFF_with_orthogroups<-RPP8_NLRs[GFF_with_orthogroups, on="ID"]
GFF_with_orthogroups$TIR_domain[is.na(GFF_with_orthogroups$TIR_domain)]<-0
GFF_with_orthogroups$CC_domain[is.na(GFF_with_orthogroups$CC_domain)]<-0
GFF_with_orthogroups$RPP8_domain[is.na(GFF_with_orthogroups$RPP8_domain)]<-0

rm(all_full_BLAST)
rm(CC_NLRs)
rm(RPP8_NLRs)
rm(TIR_NLRs)
rm(HMMER_orthogroups)


fwrite(GFF_with_orthogroups, "GFF_with_orthogroups_all_taxa_20250512.csv")

####### CALCULATE VARIABILITY IN PENNYCRESS NLRS #######

##### Import peptide seqs and select NLRs #####

# read in peptides
fs<-list.files(path="DataOutputs/peptide", pattern = ".fa", full.names = T)
fs<-fs[2:8]
pepsw <- lapply(fs, readAAStringSet)
names(pepsw) <- NULL #need to make sure it's unnamed 
ssout <- BiocGenerics::do.call(c, pepsw) #regular do.call does not let you concatenate string sets, no idea why 

# read in combBed and NLR GFF and select OGs of interest

combBed<-fread("DataOutputs/combBed.txt")
#NLR_GFF<-fread("/Users/joannarifkin/Desktop/Pennycress/NLRs/all_genes_with_positions_orthogroups_domains_HV_20241114.csv")
#fix this so it's only grabbing the NLRs 
combBedNLRs<-combBed[id %in% GFF_with_orthogroups$ID] #Subset combbed file to only retain filtered NLRs 
spl <- split(subset(combBedNLRs[genome!="TAIR11"], globOG %in% unique(GFF_with_orthogroups[genome!="TAIR11"]$Orthogroup)), by = "globOG") #split combBed file by orthogroup - datatable for every orthogroup of all the genes in that orthogroup. Select pennycress genes only by excluding Arabidopsis
# build input fastas by orthogroup


ogMd <- data.table(og = unique(GFF_with_orthogroups[genome!="TAIR11"]$Orthogroup)) #metadata file makes a fasta for every OG of interest 
fastawd="DataOutputs/orthogroup_fastas"
alnwd="DataOutputs/orthogroup_alignments"
ogMd[,`:=`(faFile = file.path(fastawd, sprintf("%s_cat.fa", og)), alFile = file.path(alnwd, sprintf("%s_aln.fa", og)))] #metadata to get file for input for mafft
ogMd[,`:=`(mafftCom = sprintf("%s > %s", faFile, alFile))]  #matrix of unaligned file and aligned files 

# create fasta files and mafft alignments for every pennycress orthogroup
for(i in 1:nrow(ogMd)){
  print(i)
  print(ogMd$og[i])
  y <- ogMd[i,]
  x <- spl[[y$og]]
  ss <- ssout[x$id] #subset amino acid string set to just the genes in that og
  writeXStringSet(ss, filepath = y$faFile)
  system2("mafft", y$mafftCom)
}


# calculate per-site shannon entropy for every orthogroup

##### Calculate entropy of the clade alignments #####
#from https://github.com/krasileva-group/hvNLR/blob/master/Atha_NLRome_CladeAnalysis.R
hvSiteEntCutoff <- 1.5
MinGapFraction <- 0.1
MinGapBlockWidth <- 1
Alph_21 <- c("A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V","-")
Alph_20 <- c("A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V")
files<-list.files(path = "DataOutputs/orthogroup_alignments", full.names = T)
EntropyNG <- vector("list",length(files))
Entropy <- vector("list",length(files))
stats<-vector()

## Entropy calculation
for (i in seq_along(files)) {
  ## Read alignment
  maa <- readAAMultipleAlignment(files[[i]])
  ## Extracting folder name
  tmp<-strsplit(files[[i]], "/")[[1]][[3]]
  og_name<-(strsplit(tmp, "_")[[1]][[1]])
  print(i)
  print(paste("Looking at clade",og_name,"..."))
  ## Filter out gappy columns
  if ("-" %in% rownames(consensusMatrix(maa))){
    autoMasked <- maskGaps(maa, min.fraction = MinGapFraction, min.block.width = MinGapBlockWidth) ##KEY FILTERING PARAMETERS
    MinAli <- as(autoMasked, "AAStringSet")
  }else{MinAli<-as(maa, "AAStringSet")}
  MinAli
  ## Calculating Consensus Matrix
  (Tidy_CM<-as_tibble(t(consensusMatrix(MinAli, baseOnly = T))))
  nrow(Tidy_CM)
  ## Compensating for consensus matrix not keeping full alphabet in output
  for (a in setdiff(Alph_21,colnames(Tidy_CM))){
    vec <- as_tibble(0*(1:nrow(Tidy_CM)))
    colnames(vec) <- paste(a)
    Tidy_CM <- as_tibble(cbind(Tidy_CM,vec))
  } 
  ##Selecting relevant columns
  (Tidy_CM_Gaps <- select(Tidy_CM,Alph_21))
  (Tidy_CM_NoGaps <- select(Tidy_CM,Alph_20))
  
  ##Entropy Calculation
  ent <- apply(Tidy_CM_Gaps, 1, entropy,unit="log2") %>% as_tibble()
  colnames(ent)<-paste0("Entropy_",og_name)
  
  
  ##Entropy Calculation Ignoring Gaps
  entNG <- apply(Tidy_CM_NoGaps, 1, entropy,unit="log2") %>% as_tibble()
  colnames(entNG)<-paste0("EntropyNoGaps_",og_name)
  #View(entNG)
  
  ##Save fraction of invariant positions with/without gaps and number of highly variable positions (without gaps)
  frbl <- length(which(ent == 0))/nrow(ent)
  frblng <- length(which(entNG == 0))/nrow(entNG)
  nHVsites <- length(which(entNG > hvSiteEntCutoff))                          ####KEY CUTOFF PARAMETER
  stats <- rbind(stats,c(og_name,nrow(maa),frbl,frblng,nrow(ent), nHVsites))  
  
  ##Save results of entropy calculation
  Entropy[[i]] <- ent
  fwrite(ent, file = paste("DataOutputs/orthogroup_entropy/", og_name, "entropy.txt", sep=""))
  EntropyNG[[i]] <- entNG
  fwrite(ent, file = paste("DataOutputs/orthogroup_entropy/", og_name, "entropy_ng.txt", sep=""))
}


stats<-as.data.table(stats)
setnames(stats, new=c("Orthogroup", "n_genes", "fract_invar", "fract_invar_no_gaps", "n_peptides", "n_HV_peptides"))
stats$n_genes<-as.numeric(stats$n_genes)
stats$fract_invar<-as.numeric(stats$fract_invar)
stats$fract_invar_no_gaps<-as.numeric(stats$fract_invar_no_gaps)
stats$n_peptides<-as.numeric(stats$n_peptides)
stats$n_HV_peptides<-as.numeric(stats$n_HV_peptides)
stats[,HVNLR:=n_HV_peptides>=10]

rm(autoMasked); rm(ent);rm(entNG);rm(maa);rm(MinAli);rm(Tidy_CM);rm(Tidy_CM_Gaps);rm(Tidy_CM_NoGaps);rm(vec)


fwrite(stats, "DataOutputs/OG_HVNLR_stats_20250501.csv")

#stats<-fread("DataOutputs/OG_HVNLR_stats_20250501.csv")


GFF_with_orthogroups<-stats[GFF_with_orthogroups, on = "Orthogroup"]



#HMMER_orthogroups_counts<-fread("../HMMER_orthogroup_counts_2024_11_14.csv")
#View(NLR_GFF)
#NLR_GFF<-NLR_GFF[stats, on = "Orthogroup"]
HMMER_orthogroups_counts<-HMMER_orthogroups_counts[stats, on="Orthogroup"]


####### IDENTIFY NLR CLUSTERS AND SYNTENIC BLOCKS #######

#simply define NLR clusters in the pennycress genomes and then see if they inhabit syntenic blocks containing Arabidopsis NLR clusters


##### Define clusters #####
GFF_with_orthogroups<-GFF_with_orthogroups[order(file_name,seqid,start,end)]
GFF_with_orthogroups[,nearest_upstream:=start-data.table::shift(end, n = 1), by = c("file_name", "seqid")]
GFF_with_orthogroups[,nearest_downstream:=end-data.table::shift(start, n = -1), by = c("file_name", "seqid")]
GFF_with_orthogroups[,in_cluster:=FALSE]
GFF_with_orthogroups[,in_cluster:=(abs(nearest_upstream)<=50000 | abs(nearest_downstream)<=50000)]
GFF_with_orthogroups$in_cluster[is.na(GFF_with_orthogroups$in_cluster)] <- FALSE
GFF_with_orthogroups[,start_cluster:=(in_cluster==TRUE & (is.na(nearest_upstream)==TRUE | abs(nearest_upstream)>50000))]
GFF_with_orthogroups[,end_cluster:=(in_cluster==TRUE & (is.na(nearest_downstream)==TRUE | abs(nearest_downstream)>50000))]
GFF_with_orthogroups[,df := c(0, start[-1] - end[1:(.N-1)]), by = c("seqid", "file_name")] # Dstance between left and right 
GFF_with_orthogroups[,sameClstr := df <= 50e3, by = c("seqid", "file_name")] #is the one to the left the same block
GFF_with_orthogroups[,clstr := cumsum(!sameClstr), by = c("file_name", "seqid")]
GFF_with_orthogroups[ ,clstr_size := .N, by = .(file_name, seqid, clstr)]
GFF_with_orthogroups[clstr_size>1]
GFF_with_orthogroups[,NLRClusterID:=paste(file_name, seqid, clstr, sep="_")]

fwrite( (GFF_with_orthogroups[clstr_size>1][, .(howmanyblocks = length(unique(clstr))), by = c("file_name", "seqid")]),"DataOutputs/NLR_block_counts_by_chromosome_and_accession20250509.csv") 
fwrite(GFF_with_orthogroups,"DataOutputs/NLR_block_genes_20250516.csv") 

##### Identify syntenic blocks containing Arabidopsis NLR clusters #####
#read in syntenic blocks
synBlk<-fread("DataOutputs/syntenicBlock_coordinates.csv")
synBlk$width<-synBlk$endBp1-synBlk$startBp1
#pull arabidopsis clusters from NLR gff
AraClusters<-GFF_with_orthogroups[clstr_size>1 & file_name=="TAIR" & (start_cluster == T | end_cluster ==T) ]

#convert arabidopsis clusters to range
AraClustersRanges<-makeGRangesFromDataFrame(AraClusters)

#pull range of arabidopsis syntenic blocks with arabidopsis as genome 1
syntenic_ranges<-makeGRangesFromDataFrame(synBlk[genome1=="TAIR11"], keep.extra.columns = T, seqnames.field = "chr1", start.field = "startBp1", end.field = "endBp1")

#select a subset of those syntenic blocks that overlap with arabidopsis nlr clusters
syntenic_ranges_w_clusters<-subsetByOverlaps(syntenic_ranges, AraClustersRanges) # subsetByOverlaps extracts the elements in the query that overlap at least one element in the subject. These are syntenic blocks that contain an Arabidopsis NLR cluster.

#turn that back into a data table so it can be refactored into ranges in the other direction - this is now a list of syntenic blocks with arabidopsis as genome 1 that overlap arabidopsis NLR clusters (about 1/9 of total clusters)
synBlk_w_cluster<-as.data.table(syntenic_ranges_w_clusters)

fwrite(synBlk_w_cluster, "DataOutputs/syntenic_blocks_that_contain_TAIR_NLR_clusters20250501.csv")

# how many total syntenic blocks per genome vs Arabidopsis? How many with NLR clusters in Arabidopsis? 
for (genome in unique(synBlk$genome1)) {
  print(genome)
  print(nrow(synBlk[genome1=="TAIR11" & genome2==genome]))
  print(nrow(synBlk_w_cluster[genome1=="TAIR11" & genome2==genome]))
}
#raw TAIR:TAIR results include chloroplast and mitochondrial genomes 

#change the names for downstream stuff - start/end will break some of the "range" stuff
setnames(synBlk_w_cluster, c("seqnames", "start","end"), c("TAIR_chr","TAIR_start","TAIR_end"))

##### Quantify block overlap by genome #####


all_cluster_ranges<-data.table()


for (gen in unique(GFF_with_orthogroups$file_name)) {
  print(gen)
  genomeClusters<-GFF_with_orthogroups[clstr_size>1 & file_name==gen & (start_cluster == T | end_cluster ==T) ] #select NLRs in that genome only and only the first and last gene in a cluster
  genomeClusters[,clstrStart:=min(start),.(seqid,clstr)] #define start of cluster as the start of the first gene in the cluster
  genomeClusters[,clstrEnd:=max(end),.(seqid,clstr)] #define end of cluster as end of last gene in cluster
  setnames(genomeClusters, c("start","end"), c("gene_start","gene_end")) #change names to avoid GenomicRanges issues
  genomeClustersRanges<-makeGRangesFromDataFrame(genomeClusters, keep.extra.columns = T, start.field = "clstrStart", end.field = "clstrEnd") #convert to GRanges object
  genomeClustersRanges<-GenomicRanges::reduce(genomeClustersRanges,ignore.strand = TRUE) #collapse first gene / last gene duplicates, and ignore strand
  clusters_only_ranges_dt<-as.data.table(genomeClustersRanges)
  setnames(clusters_only_ranges_dt, c("seqnames", "start","end"), c("seqid", "clstrStart","clstrEnd"))
  clusters_only_ranges_dt$gen<-gen
  tmp<-unique(genomeClusters[,c("seqid", "clstrStart","clstrEnd", "clstr","clstr_size")])
  clusters_only_ranges_dt<-tmp[clusters_only_ranges_dt, on=c("seqid","clstrStart", "clstrEnd")]
  all_cluster_ranges<-rbind(all_cluster_ranges, clusters_only_ranges_dt)
  genome_synBlock<-synBlk_w_cluster[genome1=="TAIR11" & grepl(pattern = gen,x = genome2)==T] #pull syntenic blocks containing NLR clusters between focal genome and Arabidopsis 
  genome_synBlock[,switch:=startBp2>endBp2] #reorient if the start of the gene is upstream in the syntenic block
  genome_synBlock$startBp2_tmp<-NA_integer_ #prepare empty columns
  genome_synBlock$endBp2_tmp<-NA_integer_
  genome_synBlock[switch==T]$startBp2_tmp<-genome_synBlock[switch==T]$endBp2 #reverse start and end if block reversed
  genome_synBlock[switch==F]$startBp2_tmp<-genome_synBlock[switch==F]$startBp2
  genome_synBlock[switch==T]$endBp2_tmp<-genome_synBlock[switch==T]$startBp2
  genome_synBlock[switch==F]$endBp2_tmp<-genome_synBlock[switch==F]$endBp2
  genome_syntenic_ranges<-makeGRangesFromDataFrame(genome_synBlock, keep.extra.columns = T, seqnames.field = "chr2", start.field = "startBp2_tmp", end.field = "endBp2_tmp") #These are only blocks syntenic between TAIR11 and the focal genome, positioned by focal genome
  genomeSyntenicClustersRanges<-subsetByOverlaps(genomeClustersRanges, genome_syntenic_ranges) # subsetByOverlaps extracts the elements in the query - the ranges containing NLR clusters in the pennycress genome - that overlap at least one element in the subject - syntenic blocks containing an arabidopsis cluster. These are syntenic blocks that contain an Arabidopsis NLR cluster.
  
  #summarise output
  print(paste( length(genomeSyntenicClustersRanges@ranges), "clusters overlapping out of", length(genomeClustersRanges@ranges), "total clusters"))
  print(paste(sum(genomeSyntenicClustersRanges@ranges@width), "base pairs in NLR clusters"))
  print("cluster size quantile")
  print(quantile(genomeSyntenicClustersRanges@ranges@width))
  print(quantile(genomeClusters$clstr_size))
}

fwrite(all_cluster_ranges, "DataOutputs/positions_of_all_NLR_clusters20250509.csv")

##### Pull Arabidopsis syntenic hits for NLR genes #####
#GFF_with_orthogroups<-fread("DataOutputs/NLR_block_genes_20250430.csv") 

synHitsfiles<-list.files(path = "DataOutputs/synhits/", pattern = "vs_TAIR11.synHits.txt.gz", full.names = T)

names(synHitsfiles)<-substr(synHitsfiles,22,40)
names(synHitsfiles)<-lapply((strsplit(names(synHitsfiles), split = "_")), `[[`, 1)
synHits<-rbindlist(lapply(names(synHitsfiles),function(i)
  data.table(genome=i,
             fread(synHitsfiles[i]))))

#synHits<-synHits[,c("genome1", "genome2", "id1","id2", "chr2", "start2", "end2","sameOG" )]
#setnames(synHits, c("id1","id2","chr2", "start2", "end2","sameOG"), c("ID","TAIR_ID","TAIR_chr", "TAIR_start", "TAIR_end","sameOG" ) )

setnames(synHits, "id1", "ID")
intersect(synHits$ID,GFF_with_orthogroups$ID)

#GFF_with_orthogroups<-synHits[GFF_with_orthogroups, on = .(ID)]
synHits_with_orthogroups<-GFF_with_orthogroups[synHits, on=.(ID, genome)]

GFFsynHits_with_orthogroups<-synHits[GFF_with_orthogroups, on=.(ID, genome)]


##### SUMMARIES AND COUNTS #####




#length(unique(GFFsynHits_with_orthogroups$ID))

#note that NLR numbers are inflated because some genes have multiple syntenic hits in Arabidopsis
for (gen in unique(synHits_with_orthogroups$genome1)) {
  print(gen) 
  GFFOGs<-rbind(table(synHits_with_orthogroups[ID %in% GFF_with_orthogroups$ID & genome1==gen]$sameOG),
              table(synHits[ID %notin% GFF_with_orthogroups$ID & genome1==gen]$sameOG))
  print(GFFOGs)
  print(chisq.test(GFFOGs))
  print(c("ratio same OG NLR", (GFFOGs[3]/(GFFOGs[3] + GFFOGs[1]))))
  print(c("ratio same OG non-NLR", (GFFOGs[4]/(GFFOGs[2] + GFFOGs[4]))))
}





GFF_with_orthogroups[,is_arabidopsis_HVNLR:=id2 %in% Arabidopsis_hvNLRs]
fwrite(GFF_with_orthogroups_and_hvNLRs, "GFF_with_orthogroups_and_hvNLRs20250508.csv")
fwrite(GFF_with_orthogroups_and_hvNLRs_and_synHits, "GFF_with_orthogroups_and_hvNLRs_and_synHits20250508.csv")

GFF_with_orthogroups_and_hvNLRs_and_synHits<-stats[GFFsynHits_with_orthogroups, on="Orthogroup"]

GFF_with_orthogroups_and_hvNLRs<-fread("DataOutputs/GFF_with_orthogroups_all_taxa_20250512.csv")
##### Summary statistics for NLRs and N-terminal types summary statistics #####

table(GFF_with_orthogroups$genome)
table(GFF_with_orthogroups$genome, GFF_with_orthogroups$TIR_domain)
table(GFF_with_orthogroups$genome, GFF_with_orthogroups$CC_domain)
table(GFF_with_orthogroups$genome, GFF_with_orthogroups$RPP8_domain)


table(GFF_with_orthogroups$includes_arabidopsis_hvNLR, GFF_with_orthogroups$RPP8_domain)
table(GFF_with_orthogroups$includes_arabidopsis_hvNLR, GFF_with_orthogroups$TIR_domain)
table(GFF_with_orthogroups$includes_arabidopsis_hvNLR, GFF_with_orthogroups$CC_domain)

table(GFF_with_orthogroups$RPP8_domain)
table(GFF_with_orthogroups$TIR_domain)
table(GFF_with_orthogroups$CC_domain)

length(unique(GFF_with_orthogroups[TIR_domain==1]$Orthogroup))
length(unique(GFF_with_orthogroups[RPP8_domain==1]$Orthogroup))
length(unique(GFF_with_orthogroups[CC_domain==1]$Orthogroup))
length(unique(GFF_with_orthogroups[CC_domain==0 & RPP8_domain==0 & TIR_domain == 0]$Orthogroup))


length(unique(GFF_with_orthogroups[genome!="TAIR11"]$Orthogroup))
length(unique(GFF_with_orthogroups[genome!="TAIR11" & TIR_domain==1]$Orthogroup))
length(unique(GFF_with_orthogroups[genome!="TAIR11" & TIR_domain==1]$Orthogroup))
length(unique(GFF_with_orthogroups[genome!="TAIR11" & RPP8_domain==1]$Orthogroup))
length(unique(GFF_with_orthogroups[genome!="TAIR11" & CC_domain==1]$Orthogroup))
length(unique(GFF_with_orthogroups[genome!="TAIR11" & CC_domain==0 & RPP8_domain==0 & TIR_domain == 0]$Orthogroup))


table(GFF_with_orthogroups$file_name, GFF_with_orthogroups$TIR_domain)
table(GFF_with_orthogroups$file_name, GFF_with_orthogroups$CC_domain)
table(GFF_with_orthogroups$file_name, GFF_with_orthogroups$RPP8_domain)


table(HMMER_orthogroups_counts$variable_size_in_pennycress)
sum(as.numeric(HMMER_orthogroups_counts$TAIR==0))
sum(as.numeric(HMMER_orthogroups_counts$variable_size_in_pennycress==1))

hist(HMMER_orthogroups_counts$size_in_pennycress)


##### Summarize and plot entropy ##### 


cor.test(stats$n_genes, stats$n_HV_peptides)
cor.test(stats[n_genes>6]$n_genes, stats[n_genes>6]$n_HV_peptides)
cor.test(stats[n_genes>9]$n_genes, stats[n_genes>9]$n_HV_peptides)
cor.test(stats[n_genes>12]$n_genes, stats[n_genes>12]$n_HV_peptides)

tmp<-HMMER_orthogroups_counts[n_genes>6]
tmp<-HMMER_orthogroups_counts[TAIR11!=0]

table(stats$HVNLR)
View(HMMER_orthogroups_counts[HVNLR==TRUE & TAIR11!=0])
table(HMMER_orthogroups_counts$HVNLR, HMMER_orthogroups_counts$TAIR11)
chisq.test(table(HMMER_orthogroups_counts$HVNLR, HMMER_orthogroups_counts$includes_arabidopsis_hvNLR))
chisq.test(table(tmp$HVNLR, tmp$includes_arabidopsis_hvNLR))

View(HMMER_orthogroups_counts[includes_arabidopsis_hvNLR==T])
View(HMMER_orthogroups_counts[HVNLR==T])

table(HMMER_orthogroups_counts$HVNLR)
mean(HMMER_orthogroups_counts[includes_arabidopsis_hvNLR==T]$n_genes)
mean(HMMER_orthogroups_counts[includes_arabidopsis_hvNLR==F]$n_genes)

hist(stats$n_HV_peptides)


plot(stats$n_genes, stats$n_HV_peptides)
plot(stats$n_peptides, stats$n_HV_peptides)

# Plot entropy for all orthogroups
pdf("FinalFigures/EntropyPlots20250512.pdf")
for (i in 1:length(Entropy)) {
  Ent <- as_tibble(cbind(rep(times=nrow(Entropy[[i]]), names(Entropy[[i]])), 1:nrow(Entropy[[i]]),Entropy[[i]]))
  colnames(Ent)<-c("OG", "Position","Entropy")
  #print(i)
  #print(Ent)
  #print(strsplit(Ent$OG[1], "_")[[1]][2] %in% HMMER_orthogroups_counts[HVNLR==T]$Orthogroup)
  hvstatus=fcase(
    strsplit(Ent$OG[1], "_")[[1]][2] %in% HMMER_orthogroups_counts[HVNLR==T]$Orthogroup, "HVNLR",
    default="")
  #print(hvstatus)
  p<-ggplot(Ent, aes(x = Position))+
    geom_line(aes(y = Entropy), color = "blue") +
    xlab("Residue") + 
    ylab("Shannon Entropy") +
    ggtitle(paste(hvstatus, "Orthogroup", strsplit(Ent$OG[1], "_")[[1]][2], "entropy")) +
    ylim(0,3) +
    theme_classic()
  plot(p)
}
dev.off()



####### SUMMARIZE CLUSTER CONTENTS #######
all_cluster_ranges<-fread("DataOutputs/positions_of_all_NLR_clusters20250509.csv") #Import clusters as ranges
GFF_with_orthogroups<-fread("DataOutputs/NLR_block_genes_20250516.csv") #Import GFF

#clusterPlot<-GFF_with_orthogroups[clstr_size>1,c("file_name", "seqid", "NLRClusterID","clstr_size", "ID", "Orthogroup","HVNLR", "start_cluster","end_cluster", "clstr", "start","end", "RPP8_domain","CC_domain","TIR_domain")]
clusterPlot<-GFF_with_orthogroups[,c("file_name", "seqid", "NLRClusterID","clstr_size", "ID", "Orthogroup","HVNLR", "includes_arabidopsis_hvNLR", "start_cluster","end_cluster", "clstr", "start","end", "RPP8_domain","CC_domain","TIR_domain")]
setnames(clusterPlot, "file_name","gen")

#Connect cluster IDs with gene IDs
tmp<-clusterPlot[,c("NLRClusterID","seqid","gen", "clstr","start")]
tmp[,clstrStart:=start]
tmp<-tmp[all_cluster_ranges[,-"strand"], on=c("seqid", "gen", "clstr" , "clstrStart")]

intersect(colnames(tmp), colnames(clusterPlot))
tmp<-tmp[,-"start"]

clusterPlot<-unique(tmp[clusterPlot, on=c("NLRClusterID", "seqid","gen", "clstr", "clstr_size"), allow.cartesian = T])
clusterPlot[is.na(clstrStart)]$clstrStart<-clusterPlot[is.na(clstrStart)]$start
clusterPlot[is.na(clstrEnd)]$clstrEnd<-clusterPlot[is.na(clstrEnd)]$end
clusterPlot[,midCluster:= (clstrStart+clstrEnd)/2]
clusterPlot[,nTerminal:=fcase(
  RPP8_domain==1, "RPP8",
  TIR_domain==1,"TIR",
  CC_domain==1,"CC",
  default = "other"
)]

fwrite(clusterPlot, "DataOutputs/OrthogroupClustersForPlotting20250516.csv")


##### Plot NLR clusters in synteny ##### 
clusterPlot<-fread("DataOutputs/OrthogroupClustersForPlotting20250512.csv")

clusterPlotOG<-clusterPlot[gen %in% c("TAIR","MN106"),c("gen", "seqid", "Orthogroup","clstr","midCluster","HVNLR", "includes_arabidopsis_hvNLR")]
setnames(clusterPlotOG, 
         c("gen", "seqid", "Orthogroup","clstr","midCluster","HVNLR"),
         c("facetGT", "facetCHR", "colorOG","factorCLSTR","positionMIDCLSTR","HVNLR"))

pdf("FinalFigures/barplotWithOrthogroups.pdf")
ggplot(transform(clusterPlotOG, 
                 gen=factor(facetGT, levels=c("AK34W"   , "Ames",     "Loretto"  ,"MN106",    "MN134"    ,"PI650286", "Tibet" , "TAIR" ))))+ aes(x = factorCLSTR, fill = colorOG) + 
  geom_bar()+
  facet_wrap(~gen+facetCHR, ncol = 7)+
  theme_light()+
  theme(legend.position = "none")
dev.off()

fwrite(clusterPlotOG, "dataForOrthogroupBarPlots20250516.csv")



clusterPlotNTermDomain<-clusterPlot[gen %in% c("TAIR","MN106"),c("gen", "seqid", "nTerminal","clstr","midCluster","HVNLR")]
setnames(clusterPlotNTermDomain, 
         c("gen", "seqid", "nTerminal","clstr","midCluster","HVNLR"),
         c("facetGT", "facetCHR", "colorNTERM","factorCLSTR","positionMIDCLSTR","HVNLR"))

fwrite(clusterPlotNTermDomain, "dataForNTerminalBarPlots20250516.csv")

pdf("FinalFigures/barplotWithNTerminals.pdf")
ggplot(transform(clusterPlotNTermDomain, 
                 gen=factor(facetGT, levels=c("AK34W"   , "Ames",     "Loretto"  ,"MN106",    "MN134"    ,"PI650286", "Tibet" , "TAIR" ))))+ aes(x = factorCLSTR, fill = colorNTERM) + 
  geom_bar()+
  facet_wrap(~facetGT+facetCHR, ncol = 7)+
  theme_light()+
  theme(legend.position = "none")
dev.off()

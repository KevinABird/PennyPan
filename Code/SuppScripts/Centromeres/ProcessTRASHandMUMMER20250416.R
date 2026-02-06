library(data.table)
library(Biostrings)
library(ggplot2)
library(stringr)
library(pwalign)

setwd("~/FinalUploads/Centromeres/")


##### Import and plot windows from TRASH output ##### 
#Import peaks file to identify main satellite sizes
peaks_files<-list.files(path = "data", pattern="*peaks.csv", recursive = T, full.names = T)
names(peaks_files)<-sapply(strsplit(peaks_files, split="\\/"),"[", 2)
peaks<-rbindlist(lapply(names(peaks_files), function(i)
  data.table(file_name=i,
             fread(file=peaks_files[i],, header = T))
))

#Import repeat summary files for all genomes
repsum_files<-list.files(path = "data", pattern="Summary", recursive = T, full.names = T)
names(repsum_files)<-sapply(strsplit(repsum_files, split="\\/"),"[", 2)
repsum<-rbindlist(lapply(names(repsum_files), function(i)
  data.table(file_name=i,
             fread(file=repsum_files[i],, header = T))
))


# Output consensus sequences
repsum<-repsum[most.freq.value.N %in% c(166,73, 523, 503, 648)] #focus on only main satellites
repsum[, file_name_match:= tstrsplit(file_name, "_", fixed=TRUE)[2]]
repsum[,center:=(start+end)/2]

consensus648<-DNAStringSet(x=(repsum[most.freq.value.N %in% c(648) & consensus.primary!="none_identified"]$consensus.primary))
names(consensus648) <-repsum[most.freq.value.N %in% c(648) & consensus.primary!="none_identified"]$consensusName
writeXStringSet(consensus648, "data/consensus648.fasta")



consensus523<-DNAStringSet(x=(repsum[most.freq.value.N %in% c(523) & consensus.primary!="none_identified"]$consensus.primary))
names(consensus523) <-repsum[most.freq.value.N %in% c(523) & consensus.primary!="none_identified"]$consensusName
writeXStringSet(consensus523, "data/consensus523.fasta")


consensus503<-DNAStringSet(x=(repsum[most.freq.value.N %in% c(503) & consensus.primary!="none_identified"]$consensus.primary))
names(consensus503) <-repsum[most.freq.value.N %in% c(503) & consensus.primary!="none_identified"]$consensusName
writeXStringSet(consensus503, "data/consensus503.fasta")


consensus166<-DNAStringSet(x=(repsum[most.freq.value.N %in% c(166) & consensus.primary!="none_identified"]$consensus.primary))
names(consensus166) <-repsum[most.freq.value.N %in% c(166) & consensus.primary!="none_identified"]$consensusName
writeXStringSet(consensus166, "data/consensus166.fasta")


consensus73<-DNAStringSet(x=(repsum[most.freq.value.N %in% c(73) & consensus.primary!="none_identified"]$consensus.primary))
names(consensus73) <-repsum[most.freq.value.N %in% c(73) & consensus.primary!="none_identified"]$consensusName
writeXStringSet(consensus73, "data/consensus73.fasta")

consensusAll<-DNAStringSet(x=c(consensusString(consensus73, ambiguityMap="N", threshold=0.5),
                               consensusString(consensus166, ambiguityMap="N", threshold=0.5),
                               consensusString(consensus503, ambiguityMap="N", threshold=0.5),
                               consensusString(consensus523, ambiguityMap="N", threshold=0.5),
                               consensusString(consensus648, ambiguityMap="N", threshold=0.5)))
names(consensusAll)<-c("consensus73","consensus166","consensus503","consensus523","consensus648")
writeXStringSet(consensusAll, "data/consensusAll.fasta")

#Plot all satellite peaks for all satellites + genomes + chromosomes
pdf("figures/TRASH_repeats_by_score_20250417.pdf", 20, 15)
p<-ggplot(repsum, aes(x=center, y=ave.score, color=as.factor(most.freq.value.N))) + 
  geom_point(aes(color=as.factor(most.freq.value.N)), alpha=0.5) + 
  theme_light() +
  facet_wrap(~ file_name_match + name, axes="all_x") +
  labs(colour = "Satellite length", x="position (Mb)", y="Proportion 1kb window")
p
dev.off()




##### Identify windows from TRASH output ##### 
#Identify regions with multiple satellite blocks within 100kb of each other
repsum<-repsum[most.freq.value.N %in% c(166,73)] #focus on only main satellites


#Output fastas from primary and secondary sequence

ConsensusPrimary<-DNAStringSet(x=(repsum[consensus.primary!="none_identified"]$consensus.primary))
names(ConsensusPrimary) <-paste(repsum[consensus.primary!="none_identified"]$file_name_match,repsum[consensus.primary!="none_identified"]$start,repsum[consensus.primary!="none_identified"]$end, sep="_")
writeXStringSet(ConsensusPrimary, "data/AllConsensusPrimary.fasta")

ConsensusSecondary<-DNAStringSet(x=(repsum[consensus.secondary!="none_identified"]$consensus.secondary))
names(ConsensusSecondary) <-paste(repsum[consensus.secondary!="none_identified"]$file_name_match,repsum[consensus.secondary!="none_identified"]$start,repsum[consensus.secondary!="none_identified"]$end, sep="_")
writeXStringSet(ConsensusSecondary, "data/AllConsensusSecondary.fasta")


ConsensusPrimary73<-DNAStringSet(x=(repsum[consensus.primary!="none_identified" & most.freq.value.N == 73]$consensus.primary))
names(ConsensusPrimary73) <-paste(repsum[consensus.primary!="none_identified"& most.freq.value.N == 73]$file_name_match,repsum[consensus.primary!="none_identified"& most.freq.value.N == 73]$start,repsum[consensus.primary!="none_identified"& most.freq.value.N == 73]$end, sep="_")
writeXStringSet(ConsensusPrimary73, "data/AllConsensusPrimary73.fasta")


ConsensusPrimary166<-DNAStringSet(x=(repsum[consensus.primary!="none_identified" & most.freq.value.N == 166]$consensus.primary))
names(ConsensusPrimary166) <-paste(repsum[consensus.primary!="none_identified"& most.freq.value.N == 166]$file_name_match,repsum[consensus.primary!="none_identified"& most.freq.value.N == 166]$start,repsum[consensus.primary!="none_identified"& most.freq.value.N == 166]$end, sep="_")
writeXStringSet(ConsensusPrimary166, "data/AllConsensusPrimary166.fasta")



View(repsum)
writeXStringSet(repsum[consensus.primary], filepath = "consensusPrimarySatellites.fasta", )  




setkey(repsum, "file_name", "name", "start","most.freq.value.N")
repsum[,nearest_upstream:=start-data.table::shift(end, n = 1), by = c("file_name", "name","most.freq.value.N")]
repsum[,nearest_downstream:=end-data.table::shift(start, n = -1), by = c("file_name", "name","most.freq.value.N")]

repsum[,in_cluster:=FALSE]
repsum[,in_cluster:=(abs(nearest_upstream)<=100000 | abs(nearest_downstream)<=100000)]
repsum[,start_cluster:=(in_cluster==TRUE & (is.na(nearest_upstream)==TRUE | abs(nearest_upstream)>=100000))]
repsum[,end_cluster:=(in_cluster==TRUE & (is.na(nearest_downstream)==TRUE | abs(nearest_downstream)>=100000))]
repsum[,df := c(0, start[-1] - end[1:(.N-1)]), by = c("file_name", "name","most.freq.value.N")] # Dstance between left and right 
repsum[,sameBlk := df <= 100e3, by = c("file_name", "name","most.freq.value.N")] #is the one to the left the same block
repsum[,blk := cumsum(!sameBlk), by = c("file_name", "name","most.freq.value.N")]
repsum[,blk_size := .N, by = .(file_name, name, blk, most.freq.value.N)]

#Identify largest cluster of satellite blocks
repsum[, max_blk_size := max(blk_size), by = c("file_name", "name")]
repsum[, is_biggest:= blk_size==max_blk_size, by= c("file_name", "name")]
repsum[, blk_score:= mean(ave.score), by= c("file_name", "name", "most.freq.value.N", "blk")]
repsum[, blk_sum_consensus_count:= sum(consensus.count), by= c("file_name", "name", "most.freq.value.N", "blk")]

tmp_table<-data.table(matrix(NA_integer_, nrow=0, ncol = 7))
setnames(tmp_table, c("filename_biggest","chr_biggest", "satellite_biggest", "block_size_biggest", "score_biggest", "start_biggest", "end_biggest"))

#Loop through satellite blocks to extract the biggest block 
for (file in unique(repsum$file_name_match)) {
  for (chr in unique(repsum$name)) {
  print(file)
  print(chr)
  tmpval_start<-repsum[is_biggest==T & start_cluster==T & file_name_match==file & name==chr]$start
  start_biggest<-rep(tmpval_start, nrow(repsum[file_name_match==file & name==chr]))
  tmpval_end<-repsum[is_biggest==T & end_cluster==T & file_name_match==file & name==chr]$end
  satellite_in_biggest<-repsum[is_biggest==T & end_cluster==T & file_name_match==file & name==chr]$most.freq.value.N
  block_size<-repsum[is_biggest==T & end_cluster==T & file_name_match==file & name==chr]$blk_size
  score<-mean(repsum[is_biggest==T & file_name_match==file & name==chr]$ave.score)
  block_size_biggest<-rep(block_size, nrow(repsum[file_name_match==file & name==chr]))
  score_biggest<-rep(score, nrow(repsum[file_name_match==file & name==chr]))  
  end_biggest<-rep(tmpval_end, nrow(repsum[file_name_match==file & name==chr]))
  chr_biggest<-rep(chr, nrow(repsum[file_name_match==file & name==chr]))
  filename_biggest<-rep(file, nrow(repsum[file_name_match==file & name==chr]))
  satellite_biggest<-rep(satellite_in_biggest, nrow(repsum[file_name_match==file & name==chr]))
  tmp_joined<-cbind(filename_biggest, chr_biggest, satellite_biggest, block_size_biggest, score_biggest, start_biggest, end_biggest)
  tmp_table<-rbind(tmp_table, tmp_joined)
  }
}

repsum$start_biggest<-as.numeric(tmp_table$start_biggest)
repsum$end_biggest<-as.numeric(tmp_table$end_biggest)

repsum[,start_region:=start_biggest-10000000]
repsum[,end_region:=end_biggest+10000000]
colnames(repsum)

fwrite(repsum, file="all_processed_repeats_20250417.csv")

candidate_regions<-unique(repsum[is_biggest==T, c("name", "start_biggest","end_biggest", "file_name_match",  "most.freq.value.N", "blk", "blk_size", "blk_score", "blk_sum_consensus_count", "is_biggest")])

# Manually relocate candidate region for MN106
tmp<-repsum[file_name_match=="MN106" & name=="Chr06" & blk_size == 11]
tmp$start_biggest<-tmp[start_cluster==T]$start
tmp$end_biggest<-tmp[end_cluster==T]$end
tmp_row<-unique(tmp[, c("name", "start_biggest","end_biggest", "file_name_match",  "most.freq.value.N", "blk", "blk_size", "blk_score", "blk_sum_consensus_count", "is_biggest")])
candidate_regions[file_name_match=="MN106" & name=="Chr06"]<-tmp_row
setkey(candidate_regions, "file_name_match", "name")
fwrite(candidate_regions, "data/Thlaspi_centromere_candidate_regions20250416.csv")



#Plot satellite peaks for all genomes + chromosomes
pdf("figures/TRASH_repeats_by_score_with_TRASH_region20250417.pdf", 20, 15)

p<-ggplot(repsum, aes(x=center, y=ave.score, color=as.factor(most.freq.value.N))) + 
  geom_point(aes(color=as.factor(most.freq.value.N)), alpha=0.5) + 
  theme_light() + 
  geom_segment(aes(x=start_biggest, y=10, xend=end_biggest, yend=10), inherit.aes = F) +
  facet_wrap(~ file_name_match + name, axes="all_x")+
  labs(colour = "Satellite length", x="position (Mb)", y="Proportion 1kb window", title = "Satellite positions with TRASH regions")
p
dev.off()

pdf("figures/TRASH_repeats_by_score_with_expandedTRASH_region20250417.pdf", 20, 15)

p<-ggplot(repsum, aes(x=center, y=ave.score, color=as.factor(most.freq.value.N))) + 
  geom_point(aes(color=as.factor(most.freq.value.N)), alpha=0.5) + 
  theme_light() + 
  geom_segment(aes(x=start_region, y=10, xend=end_region, yend=10), inherit.aes = F) + 
  geom_segment(aes(x=start_biggest, y=10, xend=end_biggest, yend=10), color="red", inherit.aes = F) +
  facet_wrap(~ file_name_match + name, axes="all_x")+
  labs(colour = "Satellite length", x="position (Mb)", y="Proportion 1kb window", title = "Satellite positions with expanded TRASH regions")
p
dev.off()

#Plot CentIER regions 

#Import CentIER centromere ranges
centier_range_files<-list.files(path = "data", pattern="*centromere_range.txt", recursive = T, full.names = T)
names(centier_range_files)<-sapply(strsplit(centier_range_files, split="\\/"),"[", 2)
centier_ranges<-rbindlist(lapply(names(centier_range_files), function(i)
  data.table(file_name=i,
             fread(file=centier_range_files[i],, header = F))
))
centier_ranges[, V4:= tstrsplit(file_name, "_", fixed=TRUE)[2]]
centier_ranges<-centier_ranges[,2:5]
setnames(centier_ranges, old=c("V1", "V2", "V3", "V4"), new=c("name","cent_start","cent_end","file_name_match"))

setkey(repsum, "file_name_match", "name")
setkey(centier_ranges, "file_name_match", "name")
repsum<-centier_ranges[repsum]

pdf("figures/TRASH_repeats_by_score_with_CentIER_region20250417.pdf", 20, 15)

p<-ggplot(repsum, aes(x=center, y=ave.score, color=as.factor(most.freq.value.N))) + 
  geom_point(aes(color=as.factor(most.freq.value.N)), alpha=0.5) + 
  theme_light() + 
  geom_segment(aes(x=cent_start, y=10, xend=cent_end, yend=10), inherit.aes = F) + 
  facet_wrap(~ file_name_match + name, axes="all_x")+
  labs(colour = "Satellite length", x="position (Mb)", y="Proportion 1kb window", title = "Satellite positions with CentIER regions")
p
dev.off()


##### Extract relevant sequence sections for MUMMER ######



# make fastas from TRASH positions 

candidate_regions

genomes<-list.files(path=Sys.glob("~Pennycress/*/*/assembly"), recursive = T, pattern=".fa.gz", full.names = T)
genomes<-genomes[2:8]
names(genomes)<-c("AK34W", "Ames32873", "LorettoMN", "MN106", "MN134", "PI650286", "Tibet33")
candidate_regions[,start_region:=start_biggest-10000000]
candidate_regions[,end_region:=end_biggest+10000000]

View(candidate_regions)  
i=1
j=1
for (i in 1:length(genomes)){
  print(genomes[i])
  genome<-readDNAStringSet(genomes[i])
  print(names(genome))
  reg<-candidate_regions[file_name_match==names(genomes)[i]]
  #subset_repsum<-repsum[file_name_match==names(genomes)[i]]
  for (j in 1:7){
    print(j)
    cenregion <- subseq(genome[j], 
                start=(reg[name==(unique(reg$name)[j])]$start_region)[1], 
                end=(reg[name==(unique(reg$name)[j])]$end_region)[1])
    print(cenregion)
    writeXStringSet(cenregion, filepath = paste("TRASH_centregion_", names(genomes)[i], substr(names(genome[j]),1,5), ".fasta", sep=""))  
  }
}


#Extract original MN106 centromere candidate region
tmp<-repsum[file_name_match=="MN106" & name=="Chr06" & blk_size == 12]
tmp$start_biggest<-tmp[start_cluster==T]$start
tmp$end_biggest<-tmp[end_cluster==T]$end
tmp_row<-unique(tmp[, c("name", "start_biggest","end_biggest", "file_name_match",  "most.freq.value.N", "blk", "blk_size", "blk_score", "blk_sum_consensus_count", "is_biggest")])
tmp_row[,start_region:=start_biggest-10000000]
tmp_row[,end_region:=end_biggest+10000000]
genome<-readDNAStringSet(genomes[4])
cenregion<-subseq(genome[6], 
       start=(tmp_row$start_region)[1], 
       end=(tmp_row$end_region)[1])
writeXStringSet(cenregion, filepath="TRASH_centregion_MN106Chr06alt.fasta")


#conda activate MUMMER

#for genome in AK34W Ames32873 LorettoMN MN106 MN134 PI650286 Tibet33
#do
#for chr in Chr01 Chr02 Chr03 Chr04 Chr05 Chr06 Chr07
#do
#nucmer --nosimplify TRASH_${genome}${chr}.fasta TRASH_${genome}${chr}.fasta -p TRASH_${genome}${chr}_self
#done
#done

##### Plot MUMMER results ##### 
#https://jmonlong.github.io/Hippocamplus/2017/09/19/mummerplots-with-ggplot2/

library(dplyr)
readDelta <- function(deltafile){
  lines = scan(deltafile, 'a', sep='\n', quiet=TRUE)
  lines = lines[-1]
  lines.l = strsplit(lines, ' ')
  lines.len = lapply(lines.l, length) %>% as.numeric
  lines.l = lines.l[lines.len != 1]
  lines.len = lines.len[lines.len != 1]
  head.pos = which(lines.len == 4)
  head.id = rep(head.pos, c(head.pos[-1], length(lines.l)+1)-head.pos)
  mat = matrix(as.numeric(unlist(lines.l[lines.len==7])), 7)
  res = as.data.frame(t(mat[1:5,]))
  colnames(res) = c('rs','re','qs','qe','error')
  res$qid = unlist(lapply(lines.l[head.id[lines.len==7]], '[', 2))
  res$rid = unlist(lapply(lines.l[head.id[lines.len==7]], '[', 1)) %>% gsub('^>', '', .)
  res$strand = ifelse(res$qe-res$qs > 0, '+', '-')
  res
}

filterMum <- function(df, minl=1000, flanks=1e4){
  coord = df %>% filter(abs(re-rs)>minl) %>% group_by(qid, rid) %>%
    summarize(qsL=min(qs)-flanks, qeL=max(qe)+flanks, rs=median(rs)) %>%
    ungroup %>% arrange(desc(rs)) %>%
    mutate(qid=factor(qid, levels=unique(qid))) %>% select(-rs)
  merge(df, coord) %>% filter(qs>qsL, qe<qeL) %>%
    mutate(qid=factor(qid, levels=levels(coord$qid))) %>% select(-qsL, -qeL)
}

#buffer by ~10% of hte plot outside of it
#always use coord_fixed() with dotplots

#Import regions from above (for labeling)
candidateRegions<-fread("data/Thlaspi_centromere_candidate_regions20250416.csv")

#Import nucmer alignment files
delta_files<-list.files(path="data/MUMMERAlignments", pattern="self.delta", full.names = T)

names(delta_files)<-sapply(strsplit(delta_files, split="\\_"),"[", 2)

mumgp_all<-rbindlist(lapply(names(delta_files), function(i)
  data.table(file_name=i,
             readDelta(delta_files[i]))
))

mumgp_all[,file_name_match:=tstrsplit(file_name, split="Chr",fixed=TRUE, keep=1)]
mumgp_all[,name:=qid]
mumgp_all<-mumgp_all[candidateRegions, on=c("file_name_match", "name")]

#setnames(tmp, c("file_name_match", "name"), c("file_name_cor", "qid"))
#mumgp_all<-mumgp_all[tmp, on=c("file_name_cor", "qid")]
mumgp_all[,rsAdj:=rs+start_region]
mumgp_all[,reAdj:=re+start_region]
mumgp_all[,qsAdj:=qs+start_region]
mumgp_all[,qeAdj:=qe+start_region]
mumgp_all$start_adjusted<-mumgp_all$start_biggest-mumgp_all$start_region
mumgp_all$end_adjusted<-mumgp_all$start_adjusted+(mumgp_all$end_biggest-mumgp_all$start_biggest)




#Plot nucmer alignments
pdf("figures/scaled_all_TRASH_region_chromosomes_20250418.pdf", 20, 20)

ggplot(mumgp_all, aes(x=rsAdj, xend=reAdj, y=qsAdj, yend=qeAdj, colour=strand)) + geom_segment(alpha=.1) +
  geom_point(alpha=.1,pch="." ) + facet_wrap(~file_name_match + name, axes="all_x") +
  theme_bw() + coord_fixed() +theme(strip.text.y=element_text(angle=180, size=5),
                                    legend.position.inside=c(.99,.01), legend.justification=c(1,0),
                                    strip.background=element_blank(),
                                    axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  xlim(0, 7e07)+ 
  geom_segment(aes(x=start_biggest, y=10, xend=end_biggest, yend=10), color="black", inherit.aes = F) +
  xlab('reference sequence') + ylab('assembly') + scale_colour_brewer(palette='Set1')+
  labs(title="Scaled dotplots of all centromeric regions")
dev.off()

pdf("figures/unscaled_all_TRASH_region_chromosomes_20250418.pdf", 20, 20)
ggplot(mumgp_all, aes(x=rs, xend=re, y=qs, yend=qe, colour=strand)) + geom_segment(alpha=.1) +
  geom_point(alpha=.1,pch="." ) + facet_wrap(~file_name_match + name, axes="all_x") +
  theme_bw() + coord_fixed() +theme(strip.text.y=element_text(angle=180, size=5),
                                    legend.position=c(.99,.01), legend.justification=c(1,0),
                                    strip.background=element_blank(),
                                    axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  geom_segment(aes(x=start_adjusted, y=10, xend=end_adjusted, yend=10), color="black", inherit.aes = F) +
  xlab('reference sequence') + ylab('assembly') + scale_colour_brewer(palette='Set1')+
  labs(title="Unscaled dotplots of all centromeric regions")
dev.off()





#Plot both versions of MN106 chr06

pdf("figures/MN106Chr06BothRegions_20250418.pdf", 6, 6)

test<-data.table(file_name_match="MN106",readDelta("data/TRASH_centregion_MN106Chr06alt_self.delta"))
#test[,file_name_match:=tstrsplit(file_name, split="Chr",fixed=TRUE, keep=1)]
test[,name:=qid]


tmp<-repsum[file_name_match=="MN106" & name=="Chr06" & blk_size == 12]
tmp$start_biggest<-tmp[start_cluster==T]$start
tmp$end_biggest<-tmp[end_cluster==T]$end
tmp_row<-unique(tmp[, c("name", "start_biggest","end_biggest", "file_name_match",  "most.freq.value.N", "blk", "blk_size", "blk_score", "blk_sum_consensus_count", "is_biggest")])
tmp_row[,start_region:=start_biggest-10000000]
tmp_row[,end_region:=end_biggest+10000000]


test<-tmp_row[test, on=c("file_name_match", "name")]
#test<-test[candidateRegions, on=c("file_name_match", "name")]
test[,rsAdj:=rs+start_region]
test[,reAdj:=re+start_region]
test[,qsAdj:=qs+start_region]
test[,qeAdj:=qe+start_region]

ggplot(test, aes(x=rsAdj, xend=reAdj, y=qsAdj, yend=qeAdj, colour=strand)) + geom_segment(alpha=.1) +
  geom_point(alpha=.1,pch="." )+theme_bw() + coord_fixed() +
  geom_segment(aes(x=start_biggest, y=10, xend=end_biggest, yend=10), color="black", inherit.aes = F) +
  xlim(0, 7e07)+
  xlab('reference sequence') + ylab('assembly') + scale_colour_brewer(palette='Set1')+
  labs(title="MN106 Chr06 largest cluster")

test<-data.table(file_name_match="MN106",readDelta("data/MUMMERAlignments/TRASH_MN106Chr06_self.delta"))
#test[,file_name_match:=tstrsplit(file_name, split="Chr",fixed=TRUE, keep=1)]
test[,name:=qid]

test<-candidateRegions[test, on=c("file_name_match", "name")]
#test<-test[candidateRegions, on=c("file_name_match", "name")]
test[,rsAdj:=rs+start_region]
test[,reAdj:=re+start_region]
test[,qsAdj:=qs+start_region]
test[,qeAdj:=qe+start_region]



ggplot(test, aes(x=rsAdj, xend=reAdj, y=qsAdj, yend=qeAdj, colour=strand)) + geom_segment(alpha=.1) +
  geom_point(alpha=.1,pch="." )+theme_bw() + coord_fixed() +
  geom_segment(aes(x=start_biggest, y=10, xend=end_biggest, yend=10), color="black", inherit.aes = F) +
  xlim(0, 7e07)+
  xlab('reference sequence') + ylab('assembly') + scale_colour_brewer(palette='Set1')+
  labs(title="MN106 Chr06 second largest cluster")

dev.off()

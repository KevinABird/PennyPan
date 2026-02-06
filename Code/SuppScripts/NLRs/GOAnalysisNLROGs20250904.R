library(topGO)
library(data.table)
library(stringr)
setwd("~/NLR/PangenomeGO")
library("org.At.tair.db")
library(GENESPACE)


##### Identify dispensable and unique genes #####
v<-tmp[["globHOG"]]
cnt_pav <- function(cb, og = "globHOG", nCores = 1, nearCore = .9, shell = .5){
  tmp <- data.table(cb)
  v <- tmp[[og]]
  tmp[,ogtmp := v]
  tmp[,size:=.N, by="ogtmp"]
#  tmp[,nGenomeTotal := .N, by = c("ogtmp", "genome")]
  tmp[,nGenome := uniqueN(genome), by = "ogtmp"]
  tmp[,propGenome := nGenome / uniqueN(tmp$genome)]
  tmp[,cls := ifelse(propGenome == 1, "core",
                     ifelse(nGenome == 1, "private",                    
                            ifelse(propGenome >= nearCore, "nearcore",
                                   ifelse(propGenome >= shell, "shell", "cloud"))))]
  outCls <- tmp[,list(
    nOG = uniqueN(ogtmp), 
    nGenes = uniqueN(paste(genome, id))), 
    by = "cls"]
  outNGs <- tmp[,list(
    nOG = uniqueN(ogtmp), 
    nGenes = uniqueN(paste(genome, id))), 
    by = "nGenome"]
  setkey(outNGs, nGenome)
  outCls[,`:=`(pavClass = factor(
    cls, levels = c("core", "nearcore", "shell", "cloud", "private")),
    cls = NULL)]
  setkey(outCls, pavClass)
  return(list(cntsByNgenomes = outNGs, cntsByPAVClasses = outCls,genes=tmp))
}
cb<-cbNoTAIR
cbWTAIR<-fread("~/NLRs/DataOutputs/combBed.txt")
cbNoTAIR<-cbWTAIR[genome!="TAIR11"]
PAV<-cnt_pav(cbNoTAIR)

genes<-PAV$genes
genes[,clsCSC:=cls]
genes[size==7 & nGenome==7]$clsCSC<-"coreSC"

cbTAIR<-cbWTAIR[genome=="TAIR11"]

OGclasses<-unique(genes[,c("globHOG", "cls")])

cbTAIR<-OGclasses[cbTAIR, on="globHOG"]

cbTAIR<-cbTAIR[!is.na(cls)]

##### GO term shell genes ######

geneList <- as.numeric(cbTAIR$cls=="shell")
names(geneList) <- cbTAIR$id

#BP MF CC
View(geneList)
GOdataShell <- new("topGOdata",
                         ontology = "BP",
                         allGenes = geneList,
                         geneSelectionFun = function(x)(x == 1),
                         annot = annFUN.org, mapping = "org.At.tair.db")


resultFisherShell <- runTest(GOdataShell, algorithm="classic", statistic="fisher")

View(GenTable(GOdataShell,fisher=resultFisherShell))


##### GO term private genes ######

geneList <- as.numeric(cbTAIR$cls=="private")
names(geneList) <- cbTAIR$id

#BP MF CC

GOdataPrivate <- new("topGOdata",
                   ontology = "BP",
                   allGenes = geneList,
                   geneSelectionFun = function(x)(x == 1),
                   annot = annFUN.org, mapping = "org.At.tair.db")


resultFisherPrivate <- runTest(GOdataPrivate, algorithm="classic", statistic="fisher")

View(GenTable(GOdataPrivate,fisher=resultFisherPrivate))

##### NLRs vs other genes ######

NLRGFF<-fread("~/NLRs/DataOutputs/GFF_with_orthogroups_all_taxa_20250512.csv")
NLRGFF<-NLRGFF[genome!="TAIR11"]

genes[,NLR:=globHOG %in% NLRGFF$hogID]
View(genes)
length(unique(genes$globOG))
ogs<-unique(genes[,c("globHOG", "cls", "clsCSC", "NLR")])

table(ogs$cls, ogs$NLR)
table(ogs$clsCSC, ogs$NLR)

#Core
#NLR
107/(107+25)

#Other
22977/(22977+1188+959+1966)
4114/(4114+22980)
#NLR genes were less likely to be included in core orthogroups (81% vs 85%)
    fisher.test(as.matrix(table(ogs$cls=="core", ogs$NLR)))
    #> fisher.test(as.matrix(table(ogs$cls=="core", ogs$NLR)))
    
    #Fisher's Exact Test for Count Data
    
    #data:  as.matrix(table(ogs$cls == "core", ogs$NLR))
    #p-value = 0.2298
    #alternative hypothesis: true odds ratio is not equal to 1
    #95 percent confidence interval:
    # 0.4896054 1.2118887
    #sample estimates:
    #odds ratio 
    # 0.7573386 
unique(ogs$cls)
  
    #more frequent in orthogroups present in 2-6 genomes (21% vs 12%), 
    
    fisher.test(as.matrix(table(ogs$cls %in% c("cloud", "shell"), ogs$NLR)))
    21/(111+21)
    3154/(23940+3154)
    #and less common in core orthogroups that are single-copy in all accessions (74% vs 91%).    

    fisher.test(as.matrix(table(ogs$clsCSC=="coreSC", ogs$NLR)))
    
    length(unique(NLRGFF$hogID))
    length(unique(NLRGFF$Orthogroup))
    76/(56+76)
    20909/(6185+20909)
    
#2-6 genomes
#NLR
12/132
#Other

#Private
#NLR

#Other

((6+5+15)/(6+5+15+110))


1-((1188+959+1966)/(1188+959+1966+22977))


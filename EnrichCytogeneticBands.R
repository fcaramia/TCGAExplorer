rm(list = ls())
library(piano)
library(gtools)
library(dplyr)
library(biomaRt)
output.dir = "~/Documents/PhD/GenderAnalysis/TCGA/Analysis/TCGAExpressionExplorerOutput/"


#Dataset
c1 = loadGSC("~/Documents/PhD/Data/MSigDB/c1.all.v6.1.symbols.gmt")

#Gene Lenghts
# gene.lenghts = read.csv("~/Documents/PhD/GenderAnalysis/TCGA/Analysis/All.genes.lengths.csv")
# all.gene.lenghts = gene.lenghts


all.g =  NULL
for (c in c1$addInfo[,1]){
  
  for (g in unlist(c1$gsc[c])){
    all.g = c(all.g,g)
  }
}

all.g = unique(all.g)

cytobands.genes = data.frame(HUGO_SYMBOL = all.g, CYTO_BAND = NA)

for (c in c1$addInfo[,1]){
  
  for (g in unlist(c1$gsc[c])){
    cytobands.genes[which(cytobands.genes$HUGO_SYMBOL==g),'CYTO_BAND'] = c
  }
}

# Get sizes from biomart
library(biomaRt)
mart.hs <- useMart("ensembl", "hsapiens_gene_ensembl")

res = getBM(attributes = c("hgnc_symbol", "entrezgene","transcript_length","band","start_position","end_position",
                           "transcript_start","transcript_end","chromosome_name"), filters = "hgnc_symbol", values = all.g, mart = mart.hs)

res %>% left_join(cytobands.genes,by = c("hgnc_symbol"="HUGO_SYMBOL")) -> res.cyto

#Do summaries



  
  
  
##RMAF DATA
RMAF = read.csv("~/Documents/PhD/GenderAnalysis/TCGA/Analysis/all.TCGA.ExMut.V2.csv", as.is=T)
RMAF$PATIENT_ID = gsub("TCGA\\-([[:alnum:]]{2})\\-([[:alnum:]]{4}).*","\\2", RMAF$SAMPLE)

#Clean low reads
RMAF$SUM = RMAF$REF_COUNT + RMAF$ALT_COUNT
RMAF %>% filter(SUM>10) -> RMAF
RMAF$RMAF = RMAF$ALT_COUNT/RMAF$SUM
RMAF$LOGIT.RMAF = logit(x = RMAF$RMAF,min = 0-0.1 , max = 1+0.1)

RMAF.genes = unique(RMAF$HUGO_SYMBOL)
mart.hs <- useMart("ensembl", "hsapiens_gene_ensembl",host="www.ensembl.org",path="/biomart/martservice")
RMAF.res = getBM(attributes = c("hgnc_symbol", "entrezgene","transcript_length","band","start_position","end_position",
                                "transcript_start","transcript_end","chromosome_name"), filters = "hgnc_symbol", values = RMAF.genes, mart = mart.hs)

RMAF.res %>% filter(!is.null(band)) %>% group_by(hgnc_symbol,band) %>%
  dplyr::summarise(transcript_length = max(transcript_length),entrez=min(entrezgene)) -> RMAF.res.summ

##Read  DNA VAF#####
DNA.VAF = read.csv("~/Documents/PhD/GenderAnalysis/TCGA/Analysis/reduced.all.TCGA.curated.mutations.csv", as.is=T)

#Clean low reads
DNA.VAF$SUM = DNA.VAF$DNA_ALT_COUNT + DNA.VAF$DNA_REF_COUNT
DNA.VAF %>% filter(SUM>10) -> DNA.VAF
DNA.VAF$DNA.VAF = DNA.VAF$DNA_ALT_COUNT/DNA.VAF$SUM
DNA.VAF$LOGIT.DNA.VAF = logit(x = DNA.VAF$DNA.VAF,min = 0-0.1 , max = 1+0.1)
DNA.VAF$PATIENT_ID = gsub("TCGA\\-([[:alnum:]]{2})\\-([[:alnum:]]{4}).*","\\2", DNA.VAF$SAMPLE_ID)
DNA.VAF$MUTATION_ID = paste(DNA.VAF$PATIENT_ID,DNA.VAF$CHROMOSOME,DNA.VAF$START_POSITION,DNA.VAF$TUMOR_SEQ_ALLELE2,sep = '.')

DNA.genes = unique(DNA.VAF$HUGO_SYMBOL)
mart.hs <- useMart("ensembl", "hsapiens_gene_ensembl",host="www.ensembl.org",path="/biomart/martservice")
DNA.res = getBM(attributes = c("hgnc_symbol", "entrezgene","transcript_length","band","start_position","end_position",
                                "transcript_start","transcript_end","chromosome_name"), filters = "hgnc_symbol", values = DNA.genes, mart = mart.hs)













###Merge DNA VAF and RMAF 

DNA.VAF %>% filter(MUTATION_ID%in%RMAF$MUTATION_ID) %>% 
  select(MUTATION_ID,DNA.VAF,LOGIT.DNA.VAF) %>%
  left_join(RMAF) -> COMB



RMAF %>% group_by(GENDER,HUGO_SYMBOL) %>% dplyr::summarise(N = n()) %>% filter(GENDER == 'male') -> RMAF.NUM.MALES
RMAF %>% group_by(GENDER,HUGO_SYMBOL) %>% dplyr::summarise(N = n()) %>% filter(GENDER == 'female') -> RMAF.NUM.FEMALES

RMAF.NUM.FEMALES$N.INV = 1 /  RMAF.NUM.FEMALES$N
RMAF.NUM.MALES$N.INV = 1/   RMAF.NUM.MALES$N

df.f = data.frame(row.names = RMAF.NUM.FEMALES$HUGO_SYMBOL,n = RMAF.NUM.FEMALES$N.INV)
df.m = data.frame(row.names = RMAF.NUM.MALES$HUGO_SYMBOL,n = RMAF.NUM.MALES$N.INV)

res.f = runGSA(geneLevelStats = df.f,geneSetStat = 'gsea', adjMethod = 'fdr',gsc=c1)
res.m = runGSA(geneLevelStats = df.m,geneSetStat = 'gsea', adjMethod = 'fdr',gsc=c1)

GSAsummaryTable(res.m)
GSAsummaryTable(res.f)

networkPlot(gsaRes = res.m,class = 'non',significance = 0.5)
networkPlot(gsaRes = res.f,class = 'non',significance = 0.5)

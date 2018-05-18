rm(list = ls())
library(dplyr)
library(ggplot2)
library(gtools)
source('ExpressionPlotsFunctions.R')
#TCGA Expression plots
#Datasets and directories and variables
datasets = c("ESCA","HNSC","LUSC","BLCA","LIHC","STAD","LGG","COAD","PAAD","READ","SKCM")
#datasets = c("READ")
output.dir = "~/Documents/PhD/GenderAnalysis/TCGA/Analysis/TCGAExpressionExplorerOutput/"
input.dir = "~/Documents/PhD/Data/TCGA_2016_01_28_BROAD/CNV.Data.Minus.Germline/"
file.sufix = '.snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.seg.txt'

##RMAF DATA
RMAF = read.csv("~/Documents/PhD/GenderAnalysis/TCGA/Analysis/all.TCGA.ExMut.csv", as.is=T)
RMAF$PATIENT_ID = gsub("TCGA\\-([[:alnum:]]{2})\\-([[:alnum:]]{4}).*","\\2", RMAF$SAMPLE)
RMAF$CNV.SEGMENT.MEAN = NA
RMAF$CNV.NPROBES = NA

for (i in datasets){
  print(i)
  #READ CNVs
  cnvs = read.csv(paste(input.dir,'/',i,'/',i,file.sufix,sep = ''),sep = '\t')
  cnvs$SAMPLE_ID = gsub("TCGA\\-([[:alnum:]]{2})\\-([[:alnum:]]{4})\\-([[:alnum:]]{3}).*",
                         "TCGA\\-\\1\\-\\2\\-\\3", cnvs$Sample)
  cnvs %>%
    filter(SAMPLE_ID%in%RMAF$SAMPLE) -> CNV.CANCER
   
  RMAF %>% 
    filter(CANCER_TYPE==i&SAMPLE%in%CNV.CANCER$SAMPLE_ID) -> RMAF.CANCER
  
  prct.num = 20
  prct = floor(nrow(RMAF.CANCER)/prct.num)
  prct.it = prct
  prct.num.it = 5
  prct.num.fix = 5
  for(mut in 1:nrow(RMAF.CANCER))
  {
    if(prct.it==mut){
      print(paste(prct.num.it,'%',sep=''))
      prct.it = prct.it + prct
      prct.num.it = prct.num.it + prct.num.fix
    }
    r = RMAF.CANCER[mut,]
    c = CNV.CANCER[which(CNV.CANCER$SAMPLE_ID==r$SAMPLE&
                           CNV.CANCER$Chromosome==r$CHROM&
                           CNV.CANCER$Start<=r$POS&
                           CNV.CANCER$End>=r$POS),]
    if(nrow(c)>1){
      print('Warning')
      print(c)
    }
    if(nrow(c)==1){
      RMAF[which(RMAF$MUTATION_ID==r$MUTATION_ID),'CNV.SEGMENT.MEAN'] = c$Segment_Mean
      RMAF[which(RMAF$MUTATION_ID==r$MUTATION_ID),'CNV.NPROBES'] = c$Num_Probes
      
    }
    
  }
  
}

write.csv(x = RMAF,file = '~/Documents/PhD/GenderAnalysis/TCGA/Analysis/all.TCGA.ExMut.CNV.csv')


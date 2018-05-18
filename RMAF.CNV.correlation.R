rm(list = ls())
library(dplyr)
library(gtools)
##RMAF DATA
RMAF = read.csv("~/Documents/PhD/GenderAnalysis/TCGA/Analysis/all.TCGA.ExMut.CNV.csv", as.is=T)
RMAF$PATIENT_ID = gsub("TCGA\\-([[:alnum:]]{2})\\-([[:alnum:]]{4}).*","\\2", RMAF$SAMPLE)

#Clean low reads
RMAF$SUM = RMAF$REF_COUNT + RMAF$ALT_COUNT
RMAF %>% filter(SUM>10) -> RMAF
RMAF$RMAF = RMAF$ALT_COUNT/RMAF$SUM
RMAF$LOGIT.RMAF = logit(x = RMAF$RMAF,min = 0-0.1 , max = 1+0.1)

signature.file =  "~/Documents/PhD/GenderAnalysis/TCGA/Analysis/TP53.lists.txt"
###Read Signatures####
sig.lists = read.delim(signature.file, sep = " ")

RMAF$RMAF.scaled = abs(RMAF$RMAF-0.5)*2

cor(RMAF$RMAF.scaled,abs(RMAF$CNV.SEGMENT.MEAN), use = 'pairwise.complete.obs')
RMAF.plot = filter(RMAF,abs(CNV.SEGMENT.MEAN)>0.2&RMAF.scaled<0.99)

cor(RMAF.plot$RMAF.scaled,abs(RMAF.plot$CNV.SEGMENT.MEAN), use = 'pairwise.complete.obs')

plot(RMAF.plot$RMAF.scaled,abs(RMAF.plot$CNV.SEGMENT.MEAN),xlab = 'CN segment mean log ratio (abs)', 
      ylab='RMAF transformed')
abline(lm(RMAF.plot$RMAF.scaled~abs(RMAF.plot$CNV.SEGMENT.MEAN)),col= 'red')
for (sl in rownames(sig.lists)){
  genes = unlist(strsplit(as.character(sig.lists[sl,'GENES']),',')) 
  signature = as.vector(sig.lists[sl,'SIGNATURE'])
  
  RMAF %>% filter(HUGO_SYMBOL%in%genes) -> RMAF.aux
  print(dim(filter(RMAF.aux,!is.na(CNV.SEGMENT.MEAN))))
  print (signature)
  cor.aux = cor(RMAF.aux$LOGIT.RMAF,RMAF.aux$CNV.SEGMENT.MEAN, use = 'pairwise.complete.obs')
  print(cor.aux)  
  
}
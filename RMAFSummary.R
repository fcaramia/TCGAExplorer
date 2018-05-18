rm(list = ls())
library(dplyr)

datasets = c("ESCA","HNSC","LUSC","BLCA","LIHC","STAD","LGG","COAD","PAAD","READ","SKCM")
input.dir = "~/Documents/PhD/GenderAnalysis/TCGA/Analysis/TCGAExpressionExplorerOutput/"
output.dir = "~/Documents/PhD/GenderAnalysis/TCGA/Analysis/TCGAExpressionExplorerOutput/"

##RMAF DATA
RMAF = read.csv("~/Documents/PhD/GenderAnalysis/TCGA/Analysis/all.TCGA.ExMut.V2.csv", as.is=T)
#Clean low reads
RMAF$SUM = RMAF$REF_COUNT + RMAF$ALT_COUNT
RMAF %>% filter(SUM>10) -> RMAF
RMAF$RMAF = RMAF$ALT_COUNT/RMAF$SUM
#RMAF$LOGIT.RMAF = logit(x = RMAF$RMAF,min = 0-0.1 , max = 1+0.1)

RMAF %>% group_by(HUGO_SYMBOL,GENDER,CANCER_TYPE) %>% dplyr::summarise(N.MUTS = n(),RMAF.MED=median(RMAF)) -> aux1

RMAF %>% group_by(HUGO_SYMBOL,CANCER_TYPE) %>% dplyr::summarise(TOTAL.MUTS = n()) -> aux.tot


aux1 %>% group_by(HUGO_SYMBOL,CANCER_TYPE)%>% dplyr::summarise(DIFF=max(RMAF.MED)-min(RMAF.MED)) -> aux.delta

aux1 %>% left_join(aux.delta)  %>% left_join(aux.tot)-> aux2

aux2$DIFF = ifelse(aux2$N.MUTS==aux2$TOTAL.MUTS,aux2$RMAF.MED,aux2$DIFF)

write.csv(aux2,"~/Documents/PhD/GenderAnalysis/TCGA/Analysis/TCGAExpressionExplorerOutput/RMAFSummary.by.Cancer.csv",row.names = F)

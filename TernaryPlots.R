rm(list = ls())
library(dplyr)
library(ggrepel)
source('ExpressionPlotsFunctions.R')
source('PropensityScoresFunctions.R')
library(gtools)
library(ggtern)
library(ggplot2)
output.dir = "~/Documents/PhD/GenderAnalysis/TCGA/Analysis/TCGAExpressionExplorerOutput/AllDataSets/"

#Datasets and directories and variables
datasets = c('ALL',"ESCA","HNSC","LUSC","BLCA","LIHC","STAD","LGG","COAD","PAAD","READ","SKCM","LUAD")
#datasets = c("READ")

#Gene sizes
gene.sizes = read.csv("~/Documents/PhD/GenderAnalysis/TCGA/Analysis/All.genes.lengths.csv")

RMAF = read.csv(paste(output.dir,'Computed.All.RMAF.csv',sep = ''))
P.Coding.X.genes = read.csv("~/Documents/PhD/GenderAnalysis/Genes.Coding.X.csv")
RMAF %>% filter(ENTREZ_GENE_ID%in%P.Coding.X.genes$entrezgene) -> RMAF

#Clean low reads
RMAF$SUM = RMAF$REF_COUNT + RMAF$ALT_COUNT

RMAF %>% mutate(RMAF = ifelse(SUM>=10,ALT_COUNT/SUM,0.0)) %>% filter(!(EXPRS.GENE=='NOT.EXPRS'&SUM<10)) %>%
  filter(!(SUM<10&Z.SCORE>-4))-> RMAF

RMAF[which(RMAF$Z.SCORE<=-4&RMAF$SUM<10),'RMAF'] = 1.0
RMAF$LOGIT.RMAF = logit(x = RMAF$RMAF,min = 0-0.1 , max = 1+0.1)
#Define Expressed Mutation RMAF
RMAF$EM.RMAF = ifelse(RMAF$RMAF>=0.75,1,0)
RMAF$SM.RMAF = ifelse(RMAF$RMAF<=0.20,1,0)
RMAF$OM.RMAF = ifelse(RMAF$SM.RMAF==0&RMAF$EM.RMAF==0,1,0)



RMAF %>% filter(EM.RMAF==1) %>% group_by(HUGO_SYMBOL,GENDER) %>% dplyr::summarise(EM.TOT=n()) -> EM.TOT 
RMAF %>% filter(SM.RMAF==1) %>% group_by(HUGO_SYMBOL,GENDER) %>% dplyr::summarise(SM.TOT=n()) -> SM.TOT 
RMAF %>% filter(OM.RMAF==1) %>% group_by(HUGO_SYMBOL,GENDER) %>% dplyr::summarise(OM.TOT=n()) -> OM.TOT 


DNA.MUTS = read.csv("~/Documents/PhD/GenderAnalysis/TCGA/Analysis/full.reduced.all.raw.TCGA.curated.mutations.csv", as.is=T)
filter.out = c("Silent",'Intron','IGR',"In_Frame_Ins" ,"In_Frame_Del", "lincRNA" )

filter(DNA.MUTS,!(VARIANT_CLASSIFICATION%in%filter.out)) -> DNA.MUTS

# ###ATRX by p53 status########
# RMAF %>% filter(HUGO_SYMBOL=='ATRX') -> ATRX.Muts
# ATRX.Muts %>% group_by(GENDER,TP53_STATUS,SM.RMAF,EM.RMAF) %>% dplyr::summarise(No=n()) -> ATRX.TABLE
# write.csv(ATRX.TABLE,"~/Documents/PhD/GenderAnalysis/TCGA/Analysis/TCGAExpressionExplorerOutput/AllDataSets/ATRX.Table.csv")
# #############################


DNA.MUTS %>% filter(CHROMOSOME=='X') %>% group_by(HUGO_SYMBOL,GENDER) %>% dplyr::summarise(DNA.TOT=n()) -> DNA.TOT

DNA.TOT %>% full_join(SM.TOT) %>% 
  full_join(EM.TOT) %>% 
  full_join(OM.TOT) %>% 
  left_join(gene.sizes,by = c("HUGO_SYMBOL" = "GENE.SYMBOL")) %>% filter(CHROMOSOME=='X')-> COMB



p53.int = read.csv("~/Documents/PhD/GenderAnalysis/TP53_interactions/Candidates.2018.03.01.csv")
p53.low = p53.int$ID



COMB$TOT = COMB$EM.TOT + COMB$SM.TOT + COMB$OM.TOT
COMB$SM.TOT.PRCT =  COMB$SM.TOT*100/COMB$TOT
COMB$EM.TOT.PRCT =  COMB$EM.TOT*100/COMB$TOT
COMB$OM.TOT.PRCT =  COMB$OM.TOT*100/COMB$TOT
COMB$TP53.INT = ifelse(COMB$HUGO_SYMBOL%in%p53.low,'yes','no')

# ggplot() + 
#   #geom_point(data = COMB, aes(x = DNA.TOT,y = SM.TOT), color = 'red', size = 1) +
#   geom_point(data = COMB, aes(x = DNA.TOT,y = EM.TOT), color = 'blue', size =1) + 
#   facet_wrap(~GENDER) + labs(x = "Number DNA Muts", y = "Number EM/SM")


ggtern(data = COMB,mapping = aes(x=EM.TOT.PRCT,y=SM.TOT.PRCT,z=OM.TOT.PRCT)) + 
  geom_point(aes(color = GENDER,shape=TP53.INT)) +
  scale_color_manual(values=c("pink1", "lightskyblue"))

rm(OM.TOT,DNA.TOT,SM.TOT,EM.TOT)
attach(COMB[which(COMB$GENDER=='female'),])
plot(SM.TOT,DNA.TOT,pch=19)
cor(DNA.TOT,EM.TOT, use = 'complete.obs')
attach(COMB[which(COMB$GENDER=='male'),])
cor(DNA.TOT,OM.TOT, use = 'complete.obs')

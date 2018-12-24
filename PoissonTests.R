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
RMAF %>% filter(HUGO_SYMBOL%in%P.Coding.X.genes$hgnc_symbol) -> RMAF

#Clean low reads
RMAF$SUM = RMAF$REF_COUNT + RMAF$ALT_COUNT

RMAF %>% mutate(RMAF = ifelse(SUM>=10,ALT_COUNT/SUM,0.0)) %>% filter(!(EXPRS.GENE=='NOT.EXPRS'&SUM<10)) %>%
  filter(!(SUM<10&Z.SCORE>-4))-> RMAF

RMAF[which(RMAF$Z.SCORE<=-4&RMAF$SUM<10),'RMAF'] = 1.0
RMAF$LOGIT.RMAF = logit(x = RMAF$RMAF,min = 0-0.1 , max = 1+0.1)
#Define Expressed Mutation RMAF
RMAF$EM.RMAF = ifelse(RMAF$RMAF>=0.75,1,0)
RMAF$SM.RMAF = ifelse(RMAF$RMAF<=0.20,1,0)

RMAF %>% filter(EM.RMAF==1) %>% group_by(HUGO_SYMBOL,GENDER) %>% dplyr::summarise(EM.TOT=n()) -> EM.TOT 
RMAF %>% filter(SM.RMAF==1) %>% group_by(HUGO_SYMBOL,GENDER) %>% dplyr::summarise(SM.TOT=n()) -> SM.TOT 

DNA.MUTS = read.csv("~/Documents/PhD/GenderAnalysis/TCGA/Analysis/full.reduced.all.raw.TCGA.curated.mutations.csv", as.is=T)
DNA.MUTS %>% filter(CHROMOSOME=='X') %>% group_by(HUGO_SYMBOL,GENDER) %>% dplyr::summarise(DNA.TOT=n()) -> DNA.TOT

DNA.TOT %>% full_join(SM.TOT) %>% 
  full_join(EM.TOT) %>% 
  left_join(gene.sizes,by = c("HUGO_SYMBOL" = "GENE.SYMBOL")) %>% filter(CHROMOSOME=='X')-> COMB

ggplot() + 
  #geom_point(data = COMB, aes(x = DNA.TOT,y = SM.TOT), color = 'red', size = 1) +
  geom_point(data = COMB, aes(x = DNA.TOT,y = EM.TOT), color = 'blue', size =1) + 
  facet_wrap(~GENDER) + labs(x = "Number DNA Muts", y = "Number EM/SM")

ggtern(data = COMB,mapping = aes(x=EM.TOT,y=SM.TOT,z=DNA.TOT)) + geom_point(aes(fill = GENDER), 
                                                                             size = 1, shape =21)


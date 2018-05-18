library(dplyr)
library(tidyr)

#Clean the workspace
rm(list = ls())
#TCGA Expression
#Read Files
datasets = c("ESCA","HNSC","LUSC","BLCA","LIHC","STAD","LGG","COAD","PAAD","READ","SKCM")
datasets = c("LIHC")
output.dir = "~/Documents/PhD/GenderAnalysis/TCGA/Analysis/TCGAExpressionExplorerOutput/"
file.name = 'RMAF.csv'
file.name.x = 'X.RMAF.csv'


#################
###Clinical Data####
demo = read.csv("~/Documents/PhD/GenderAnalysis/TCGA/Analysis/all.TCGA.curated.clinical.csv", as.is = T)
demo$PATIENT_ID = toupper(demo$PATIENT_ID)

####RMAF Data####
rmaf.data = read.csv("~/Documents/PhD/GenderAnalysis/TCGA/Analysis/all.TCGA.ExMut.csv", as.is=T)

for (i in datasets)
{
  
  dir.create(paste(output.dir,i,'RMAF',sep = "/"),showWarnings = FALSE)
  rmaf.data %>% filter(CANCER_TYPE==i) -> rmaf.cancer
 
  print(i)
  print("Nr. Samples")
  print(length(unique(rmaf.cancer$PATIENT)))
  
  rmaf.cancer$rtac = ifelse(rmaf.cancer$ALT_COUNT>rmaf.cancer$NOISE_COUNT,
                            rmaf.cancer$ALT_COUNT+rmaf.cancer$REF_COUNT,
                            rmaf.cancer$NOISE_COUNT+rmaf.cancer$REF_COUNT)
  
  #Filter mutations with at least 10 reads
  filter(rmaf.cancer,rtac>=10) -> rmaf.cancer
  
  rmaf.cancer$rmaf = ifelse(rmaf.cancer$ALT_COUNT>rmaf.cancer$NOISE_COUNT,
                            rmaf.cancer$ALT_COUNT/rmaf.cancer$rtac,
                            rmaf.cancer$NOISE_COUNT/rmaf.cancer$rtac)
  
  #If gene has multiple mutations use the one with max expression
  rmaf.cancer %>% select(HUGO_SYMBOL,SAMPLE,rmaf) %>% 
                  group_by(HUGO_SYMBOL,SAMPLE) %>% 
                  summarise(rmaf_max=max(rmaf)) -> indexes
  
  rmaf.cancer %>% 
    filter(CHROM=='X') %>%
    select(HUGO_SYMBOL,SAMPLE,rmaf) %>% 
    group_by(HUGO_SYMBOL,SAMPLE) %>% 
    summarise(rmaf_max=max(rmaf)) -> indexes.x
  
  #Use only mutations present in at least 5% the population
  
  spread(indexes,HUGO_SYMBOL,rmaf_max) -> res   
  spread(indexes.x,HUGO_SYMBOL,rmaf_max) -> res.x   
  
  
  write.csv(res,paste(output.dir,i,'RMAF',file.name,sep = "/"),row.names = F, quote = F)
  write.csv(res.x,paste(output.dir,i,'RMAF',file.name.x,sep = "/"),row.names = F, quote = F)
  
}
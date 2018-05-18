rm(list = ls())
library(dplyr)

input.dir = "~/Documents/PhD/GenderAnalysis/TCGA/Analysis/TCGAExpressionExplorerOutput/"
output.dir = "~/Documents/PhD/GenderAnalysis/TCGA/Analysis/TCGAExpressionExplorerOutput/"


##RMAF DATA
RMAF = read.csv("~/Documents/PhD/GenderAnalysis/TCGA/Analysis/all.reduced.TCGA.curated.mutations.csv", as.is=T)
RMAF$PATIENT_ID = gsub("TCGA\\-([[:alnum:]]{2})\\-([[:alnum:]]{4}).*","\\2", RMAF$SAMPLE)

#Clean low reads
RMAF$SUM = RMAF$REF_COUNT + RMAF$ALT_COUNT
RMAF %>% filter(SUM>10) -> RMAF
RMAF$RMAF = RMAF$ALT_COUNT/RMAF$SUM
RMAF$LOGIT.RMAF = logit(x = RMAF$RMAF,min = 0-0.1 , max = 1+0.1)


##Do something about noise?######

###Clinical Data####
demo = read.csv("~/Documents/PhD/GenderAnalysis/TCGA/Analysis/all.TCGA.curated.clinical.csv", as.is = T)
demo$PATIENT_ID = toupper(demo$PATIENT_ID)
###############

###Propensity Scores All####
weights.all = read.csv("~/Documents/PhD/GenderAnalysis/TCGA/Analysis/TCGAExpressionExplorerOutput/AllDataSets/PS/PS.csv")

signature.file =  "~/Documents/PhD/GenderAnalysis/TCGA/Analysis/TP53.lists.txt"
###Read Signatures####
sig.lists = read.delim(signature.file, sep = " ")

RMAF %>%
  left_join(demo[c('GENDER','PATIENT_ID')]) %>%
  left_join(weights.all) %>%
  filter(!is.na(GENDER) &!is.na(WEIGHTS)) -> RMAF

RMAF$GENDER.TXT = RMAF$GENDER
dir.create(paste(output.dir,'AllDataSets',"/MUTS",sep = ""),showWarnings = F)
for (s in rownames(sig.lists)){
  
  genes = unlist(strsplit(as.character(sig.lists[s,'GENES']),',')) 
  RMAF.S = filter(RMAF,HUGO_SYMBOL%in%genes)
  signature = as.vector(sig.lists[s,'SIGNATURE'])
  print(signature)
  
  RMAF.S$GENOME = 'hg19'
  RMAF.S$MUTAF = paste(RMAF.S$GENOME,RMAF.S$CHROMOSOME,RMAF.S$START_POSITION,
                       RMAF.S$REFERENCE_ALLELE,RMAF.S$TUMOR_SEQ_ALLELE2,sep = ',')
  RMAF.S %>% 
    select(MUTAF,GENDER) -> res

  write.table(res,paste(output.dir,'AllDataSets',"/MUTS/",signature,'.tsv',sep = ''),
            sep = '\t',col.names = F, row.names = F,quote = F)
}


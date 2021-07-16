library(data.table)
library(tidyr)
library(stringr)
library(dplyr)
germs = fread('~/Documents/PhD/FilterTCGASNPs/FilteredHetGermlineBRCAXChr.tsv')

germs = germs%>% filter(str_length(REFERENCE_ALLELE)==1,str_length(TUMOR_SEQ_ALLELE2)==1)

germs = germs%>% mutate(het = as.double(gsub(".*;HET=(.*);HOM.*","\\1",INFO))) 

germs %>% group_by(START_POSITION) %>% summarise (n=n()) %>% filter(n>=50) -> pos_lst

germs = germs %>% filter(START_POSITION%in%pos_lst$START_POSITION)

germs = germs %>% select(CHROMOSOME,START_POSITION,REFERENCE_ALLELE,TUMOR_SEQ_ALLELE2) %>% distinct()

germs = germs %>% mutate(sample_id = "TCGAGold",var="SNP")

colnames(germs) = c("chr","pos","ref","alt","sample_id",'var')

fwrite(germs,"../FilteredTCGAXchrForExMutv6TCGAGold.csv")


library(data.table)
library(tidyr)
library(dplyr)
germs = fread('~/Documents/PhD/FilterTCGASNPs/FilteredHetGermlineBRCAXChr.tsv')

somat = fread('~/Documents/PhD/GenderAnalysis/TCGA/Analysis/Mutations.For.ExMutv5.AllChromosomes.csv')

germs %>% 
  mutate(TUMOR_SEQ_ALLELE2 = strsplit(as.character(TUMOR_SEQ_ALLELE2), ",")) %>% 
  unnest(TUMOR_SEQ_ALLELE2) -> germs

germs$SAMPLE_EXP = germs$SAMPLE_ID

germs$SAMPLE_ID = gsub("(TCGA\\-[[:alnum:]]{2}\\-[[:alnum:]]{4}-[[:alnum:]]{3}).*","\\1", germs$SAMPLE_EXP)

germs$VARIANT_TYPE = ifelse(nchar(germs$TUMOR_SEQ_ALLELE2)>1,'INDEL','SNP')

write.csv(germs,"~/Documents/PhD/GenderAnalysis/TCGA/Analysis/Germline.Het.BRCA.Normals.Xchrom.For.ExMutv6.csv", row.names = F)

germs$SAMPLE_CODE = gsub("TCGA\\-([[:alnum:]]{2})\\-([[:alnum:]]{4}\\-([[:alnum:]]{2})).*","\\3", germs$SAMPLE_ID)

germs$SAMPLE_TYPE <- ifelse(as.integer(germs$SAMPLE_CODE)>=20, 'Control',
                            ifelse(as.integer(germs$SAMPLE_CODE)<10, 'Tumor', 
                                   ifelse(as.integer(germs$SAMPLE_CODE)==10, 'Blood', 'Normal Tissue')))

aux = as.vector(somat[match(germs$PATIENT_ID,somat$PATIENT_ID),'SAMPLE_ID'])
germs$SOMATIC_SAMPLE_ID = aux$SAMPLE_ID
rm(aux)
germs$GERMLINE_SAMPLE_ID = germs$SAMPLE_ID
germs$SAMPLE_ID = germs$SOMATIC_SAMPLE_ID

write.csv(germs,"~/Documents/PhD/GenderAnalysis/TCGA/Analysis/Germline.Het.BRCA.Tumors.Xchrom.For.ExMutv6.csv", row.names = F)


rm(list=ls())
library(data.table)
muts = fread("~/Documents/PhD/Data/TCGA_Mutations_PAPER/SKCM_correctedMAF_mutations.out")
colnames(muts) = toupper(colnames(muts))
muts$SAMPLE_ID = gsub("TCGA\\-([[:alnum:]]{2})\\-([[:alnum:]]{4})\\-([[:alnum:]]{3}).*","TCGA-\\1-\\2-\\3", muts$TUMOR_SAMPLE_BARCODE)

out = muts[,c("CHROMOSOME","START_POSITION","REFERENCE_ALLELE","TUMOR_SEQ_ALLELE2","VARIANT_TYPE","SAMPLE_ID")] 

fwrite(out,"~/Documents/PhD/Data/TCGA_Mutations_PAPER/SKCM_correctedMAF_mutations.out.forExMut")


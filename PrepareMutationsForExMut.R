library(dplyr)
rm(list = ls())
muts = read.csv("~/Documents/PhD/GenderAnalysis/TCGA/Analysis/reduced.all.raw.TCGA.curated.mutations.csv")


##Change GDC mutations from hg38 to hg19 
library(rtracklayer)
path = system.file(package="liftOver", "extdata", "hg38ToHg19.over.chain")
hg38Tohg19 = import.chain(path)
range=GRanges(seqnames=paste('chr',muts$CHROMOSOME,sep = ''),ranges=IRanges(start=muts$START_POSITION,end=muts$END_POSITION),MUTATION_ID = muts$MUTATION_ID)

ch = liftOver(range,hg38Tohg19)
ch = as.data.frame(ch)


#Merge with all the mutations info
ch %>% left_join(muts) -> muts

muts$START_POSITION = muts$start
muts$END_POSITION = muts$end

muts$start = NULL
muts$end = NULL
muts$seqnames = NULL
muts$group = NULL
muts$group_name = NULL
muts$width = NULL
muts$strand = NULL


##Add Firehose high confidence mutations


firehose = read.csv("~/Documents/PhD/GenderAnalysis/TCGA/Analysis/reduced.all.firehose.TCGA.curated.mutations.csv", as.is=T)
firehose$PATIENT_ID = gsub("TCGA\\-([[:alnum:]]{2})\\-([[:alnum:]]{4}).*","\\2", firehose$SAMPLE_ID)

firehose$MUTATION_ID = paste(firehose$PATIENT_ID,firehose$HUGO_SYMBOL,firehose$CHROMOSOME,firehose$START_POSITION,
                         firehose$END_POSITION,firehose$TUMOR_SEQ_ALLELE2,sep = '.')

firehose = firehose[!duplicated(firehose$MUTATION_ID),]

muts$MUTATION_ID = paste(muts$PATIENT_ID,muts$HUGO_SYMBOL,muts$CHROMOSOME,muts$START_POSITION,
                         muts$END_POSITION,muts$TUMOR_SEQ_ALLELE2,sep = '.')

firehose %>% filter(!(MUTATION_ID%in%muts$MUTATION_ID)) -> firehose.not.in.muts

firehose.not.in.muts$POLYPHEN = NA
firehose.not.in.muts$FILTER = NA
firehose.not.in.muts$CONSEQUENCE = NA
firehose.not.in.muts$SIFT = NA
firehose.not.in.muts$BIOTYPE = NA
firehose.not.in.muts$IMPACT = NA
firehose.not.in.muts$VARIANT_CALLER = "GDAC"
firehose.not.in.muts$MUTATION_ID_CALLER = paste(firehose.not.in.muts$MUTATION_ID,firehose.not.in.muts$VARIANT_CALLER,sep = '.')


demo = read.csv("~/Documents/PhD/GenderAnalysis/TCGA/Analysis/all.TCGA.curated.clinical.csv", as.is = T)
demo$PATIENT_ID = toupper(demo$PATIENT_ID)

firehose.not.in.muts %>% 
  left_join(demo[c('GENDER','PATIENT_ID')]) %>%
  filter(!is.na(GENDER)) -> firehose.not.in.muts 


rbind(muts,firehose.not.in.muts[,colnames(muts)]) -> all.muts

write.csv(all.muts,"~/Documents/PhD/GenderAnalysis/TCGA/Analysis/full.reduced.all.raw.TCGA.curated.mutations.csv", row.names = F)


datasets = c("ESCA","HNSC","LUSC","BLCA","LIHC","STAD","LGG","COAD","PAAD","READ","SKCM")
all.muts %>% filter(CHROMOSOME=='X',CANCER_TYPE%in%datasets) -> muts
write.csv(muts,"~/Documents/PhD/GenderAnalysis/TCGA/Analysis/Mutations.For.ExMutv3.csv", row.names = F)

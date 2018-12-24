rm(list = ls())
library(dplyr)
library(reshape2)
#TCGA merge mutations curations
#Read Mutations

annot = read.delim("~/Documents/PhD/Data/TCGA_Xena/illuminaMethyl450_hg19_GPL16304_TCGAlegacy")
demo = read.csv("~/Documents/PhD/GenderAnalysis/TCGA/Analysis/all.TCGA.curated.clinical.csv", as.is = T)
demo$PATIENT_ID = toupper(demo$PATIENT_ID)
x.probes = annot[which(annot$chrom=='chrX'),'X.id']
write(paste(x.probes),"~/Documents/PhD/Data/TCGA_Xena/Methylation450/X.probes.txt")

methyl.x.chrom = read.delim("~/Documents/PhD/Data/TCGA_Xena/Methylation450/Melted_PANCAN_HumanMethyl_450.Xchrom")

methyl.x.chrom$PATIENT_ID = gsub("TCGA\\-([[:alnum:]]{2})\\-([[:alnum:]]{4}).*","\\2", methyl.x.chrom$Sample)
methyl.x.chrom.disp.set = left_join(methyl.x.chrom,demo[,c('PATIENT_ID','CANCER_TYPE')])
methyl.x.chrom.disp.set = methyl.x.chrom.disp.set[which(!is.na(methyl.x.chrom.disp.set$CANCER_TYPE)),]

write.csv(methyl.x.chrom.disp.set,"~/Documents/PhD/Data/TCGA_Xena/Methylation450/Melted_PANCAN_HumanMethyl_450.Xchrom.disp.set.csv")

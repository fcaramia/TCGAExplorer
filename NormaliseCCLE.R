rm (list = ls())
library(data.table)
library(reshape)
library(ggplot2)
library(dplyr)
source("~/Documents/Work/analysis/RNA.Analysis/NormalisationFunctions.R")
ccle.rna = as.data.frame(fread("~/Documents/PhD/Data/CCLE/CCLE_DepMap_18q3_RNAseq_reads_20180718.gct"))
colnames(ccle.rna) = gsub("(.*)\\s.*",'\\1',colnames(ccle.rna))
##Collapse repeated genes
ccle.rna$Name = NULL
ccle.rna.raw <- aggregate(. ~ Description, data = ccle.rna, max)
rownames(ccle.rna.raw) = ccle.rna.raw$Description
ccle.rna.raw$Description = NULL
###Normalise####
ccle.rna.norm = as.data.frame(DoCPMNorm(raw.counts = ccle.rna.raw,pop.genes = .1))

ccle.rna.norm$Gene.ID = rownames(ccle.rna.norm)
write.csv(ccle.rna.norm,'~/Documents/PhD/Data/CCLE/CCLE.NORM.COUNTS.csv', row.names = F)

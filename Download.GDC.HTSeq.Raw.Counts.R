rm(list = ls())
library(TCGAbiolinks)
library("SummarizedExperiment")
library(data.table)
datasets = c("ESCA","HNSC","LUSC","BLCA","LIHC","STAD","LGG","COAD","PAAD","READ","SKCM",'KIRC','LUAD','BRCA','OV','PRAD')
datasets = c("BRCA")
output.dir = "~/Documents/PhD/Data/GDC2020/HTSeqCounts/"
dir.create(output.dir,showWarnings = F, recursive = T)


for (i in datasets){
  dir.create(paste(output.dir,i,sep = ""),showWarnings = F)
  query <- GDCquery(project = paste("TCGA-",i, sep = ''),
                    data.category = "Transcriptome Profiling",
                    data.type = "Gene Expression Quantification", 
                    workflow.type = "HTSeq - FPKM-UQ")
  GDCdownload(query = query, directory = '~/Documents/PhD/Data/GDC2020/')
  data <- GDCprepare(query, directory = '~/Documents/PhD/Data/GDC2020/')
  res <- data.frame(assay(data))
  res$Feature = as.vector(paste(rownames(res),data@rowRanges$external_gene_name,sep = '|'))
  fwrite(res,file =paste(output.dir,i,'/',i,'.FPKM-UQ.csv',sep = ''))  
  
}


##Delete Aux files in GDC2020

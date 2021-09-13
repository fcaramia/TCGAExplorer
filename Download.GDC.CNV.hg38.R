rm(list = ls())
library(TCGAbiolinks)
library("SummarizedExperiment")
library(data.table)
library(AnnotationHub)
library(dplyr)
datasets = c("ESCA","HNSC","LUSC","BLCA","LIHC","STAD","LGG","COAD","PAAD","READ","SKCM",'KIRC','LUAD','BRCA','OV','PRAD')
datasets = c("BRCA")
output.dir = "~/Documents/PhD/Data/GDC2020/CNV/"
dir.create(output.dir,showWarnings = F, recursive = T)

for (i in datasets){
  dir.create(paste(output.dir,i,sep = ""),showWarnings = F)
  query <- GDCquery(project = paste("TCGA-",i, sep = ''),
                    data.category = "Copy Number Variation",
                    data.type = "Masked Copy Number Segment"
                    )
                    
  GDCdownload(query = query, directory = '~/Documents/PhD/Data/GDC2020/')
  data <- GDCprepare(query, directory = '~/Documents/PhD/Data/GDC2020/')
 
  data$GDC_Aliquot=NULL
  #res$Feature = as.vector(paste(rownames(res),data@rowRanges$external_gene_name,sep = '|'))
  fwrite(data,file =paste(output.dir,i,'/',i,'.hg38.CNV.csv',sep = ''))  
  
}

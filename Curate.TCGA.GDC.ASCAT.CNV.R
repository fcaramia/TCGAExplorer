rm(list = ls())
library(dplyr)
library(data.table)
library(stringr)
library(AnnotationHub)
#TCGA merge mutations curations
#Read Mutations
datasets = c("ESCA","HNSC","LUSC","BLCA","LIHC","STAD","LGG","COAD","PAAD","READ","SKCM",'KIRC','LUAD','BRCA','OV','PRAD',"GBM")
datasets = c('BRCA')
dir = "~/Documents/PhD/Data/TCGA_2016_01_28_BROAD/CNV.ASCAT2/"
badchars <- "[\xb5]|[\n]|[,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[/]|[\\]|[.]|[_]|[ ]"
cnv = NULL

files = list.files(dir)
sample_sheet = fread("~/Documents/PhD/Data/TCGA_2016_01_28_BROAD/CNV.ASCAT2/gdc_sample_sheet.2021-07-26.tsv")
for (f in files)
{
  print(f)
  cancer.type = gsub("TCGA\\-([[:alnum:]]*)\\..*","\\1",f)
  if (!(cancer.type %in% datasets)){
    next
  }
  aux_data = read.delim(paste(dir,f,sep = ""), as.is = T, comment.char = '#')
  sample_id = gsub(".*(TCGA\\-[[:alnum:]]*\\-[[:alnum:]]*\\-01.).*","\\1",sample_sheet[which(sample_sheet$`File Name`==f),"Sample ID"])
  aux_data$sample_id= sample_id
  aux_data$cancer_type = cancer.type
  aux_data$GDC_Aliquot=NULL
  
  #Make colnames consistent
  
  
  if(is.null(cnv))
  {
    cnv = aux_data 
  }
  else
  {
    #Merge datasets
    cnv <- rbind(cnv, aux_data)
  }
  #clean up
  rm(aux_data)
  
}


write.csv(cnv,"~/Documents/PhD/GenderAnalysis/TCGA/Analysis/CNV_ASCAT.csv", row.names = F)









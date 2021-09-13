####Curate Uncompressed expression file from FireHose########
library(data.table)
#TCGA Expression
#Read Files
datasets = c("ESCA","HNSC","LUSC","BLCA","LIHC","STAD","LGG","COAD","PAAD","READ","SKCM","OV")
datasets = c("GBM")

dir = "~/Documents/PhD/Data/TCGA_2016_01_28_BROAD/Expression.Data/"

for (i in datasets)
{
  print(i)
  aux_data = data.table::fread(paste(dir,i,'/',i,".rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.data.txt",sep = ""))
  aux_data = as.data.frame(aux_data[-1,])
  rownames(aux_data) = aux_data$`Hybridization REF`
  aux_data = aux_data[,-1]
  aux_data = aux_data[,c(seq(1,ncol(aux_data),3))]
  write.csv(aux_data,paste(dir,i,"/RawExpression.csv",sep = ""))
  
  #clean up
  rm(aux_data)
  
}

rm(datasets,i,dir)

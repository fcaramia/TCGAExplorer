rm(list = ls())
library(dplyr)
#Datasets and directories and variables
datasets = c("ESCA","HNSC","LUSC","BLCA","LIHC","STAD","LGG","COAD","PAAD","READ","SKCM")
input.dir = "~/Documents/PhD/GenderAnalysis/TCGA/Analysis/TCGAExpressionExplorerOutput/"
output.dir = "~/Documents/PhD/GenderAnalysis/TCGA/Analysis/TCGAExpressionExplorerOutput/"
test.genes = c("NFE2L2")
disc.exprs.int = c(0,.25,.75,1) 

####Read Annotation file####
gene.annot = read.csv("~/Documents/PhD/GenderAnalysis/TCGA/Analysis/annotation.all.genes.csv",as.is = T)
gene.annot = gene.annot[,-1]
gene.annot = gene.annot[gene.annot$Chrm!='',]
gene.annot = gene.annot[order(gene.annot$Chrm,gene.annot$Entrez.Gene.ID),]
#########################

###Mutations####
mutations = read.csv("~/Documents/PhD/GenderAnalysis/TCGA/Analysis/reduced.all.raw.TCGA.curated.mutations.csv", as.is = T)
mutations$PATIENT_ID = gsub("TCGA\\-([[:alnum:]]{2})\\-([[:alnum:]]{4}).*","\\2", mutations$SAMPLE_ID)
###############

##Big loop
for(cancer in datasets){
  print(cancer)
  ###Create directories####
  dir.create(paste(output.dir,cancer,sep = ""),showWarnings = F)
  dir.create(paste(output.dir,cancer,"/ExpressionMutationCorrelation",sep = ""),showWarnings = F)
  #########################
  ###Read Norm counts#####
  print('Reading Norm counts...')
  norm.counts = read.csv(paste(input.dir,cancer,'/Normalisation/Norm.Log.Counts.csv',sep = ""), as.is = T)
  rownames(norm.counts) = norm.counts[,1]
  norm.counts = norm.counts[,-1]
  #########################
  pat.ids.type = gsub("TCGA\\.([[:alnum:]]{2})\\.([[:alnum:]]{4})\\.([[:alnum:]]{2}).*","\\2-\\3", colnames(norm.counts))
  pat.ids = gsub("TCGA\\.([[:alnum:]]{2})\\.([[:alnum:]]{4})\\.([[:alnum:]]{2}).*","\\2", colnames(norm.counts))
  pat.ids.norm = paste(pat.ids,'11',sep = '-')
  colnames(norm.counts) = gsub("TCGA\\.([[:alnum:]]{2})\\.([[:alnum:]]{4})\\.([[:alnum:]]{2}).*","\\2-\\3", 
                               colnames(norm.counts))
  cancer.counts = norm.counts[,-which(colnames(norm.counts)%in%pat.ids.norm)]
  colnames(cancer.counts) = gsub("([[:alnum:]]{4})\\-([[:alnum:]]{2}).*","\\1", colnames(cancer.counts))
  rownames(cancer.counts) = cancer.counts$X
  cancer.counts$X = NULL
  
  disc.exprs = t(apply(cancer.counts,1,function(x) {i1=quantile(x,probs=disc.exprs.int)[2]; i2=quantile(x,probs=disc.exprs.int)[3]; ifelse(x<=i1,'low',ifelse(x>=i2,'high','med'))}))
  mutations %>% filter(CANCER_TYPE==cancer) %>% cancer.muts
  disc.muts = data.frame(row.names = unique(cancer.muts$ENTREZ_GENE_ID),col.names = unique(cancer.muts$PATIENT_ID))
  
  
  
}

rm(list = ls())
library(dplyr)
library(plyr)
library(Matrix)
library(gridExtra)
library(gtable)
library(grid)
library(fgsea)
#library(xlsx)
#library(combinat)
library(data.table)
source("ExpressionPlotsFunctions.R")
datasets = c("ESCA","HNSC","LUSC","BLCA","LIHC","STAD","LGG","COAD","PAAD","READ","SKCM","LUAD")
input.dir = "/home/fcaramia/Documents/PhD/GenderAnalysis/TCGA/Analysis/TCGAExpressionExplorerOutput/"
output.dir = "/home/fcaramia/Documents/PhD/GenderAnalysis/TCGA/Analysis/TCGAExpressionExplorerOutput/"
signature.file =  "~/Documents/PhD/GenderAnalysis/TP53_interactions/Candidates.2018.03.01.csv"

sig.lists = fread(signature.file)

genes = sig.lists$ID
neg.genes = c('CDK16', 'HUWE1', 'NOX1', 'PPEF1', 'PSMD10', 'TAF1', 'UBE2A', 'UTP14A','CUL4B', 'DDX53', 'MAGEA2')
###Start Main Loop###
print('Main Loop')
dir.create(paste(output.dir,'GSEA',sep = ""),showWarnings = F)
gset = list(P53.NEG.INT=neg.genes)
#gset = list(P53.INTERACTOR.X = genes)
res = NULL
pdf(paste(output.dir,'GSEA/plots.NEG.FC.pdf',sep = ""))
for (cancer in datasets)
{
  print(cancer)
  ###Create directories####
  
  
  
  #########################
  files = list.files(path = paste(input.dir,cancer,'/DifferentialExpression/',sep = ''), 
                     pattern = '.csv', all.files = FALSE,
             full.names = FALSE, recursive = FALSE,
             ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
  

  for (f in files){

    if(file.exists(paste(input.dir,cancer,'/DifferentialExpression/',f,sep = ''))){
      
      diffexp = read.csv(paste(input.dir,cancer,'/DifferentialExpression/',f,sep = ''), as.is = T)
      m = gsub("(.*)\\.csv",'\\1',f)
      
      
      diffexp %>% 
        select_('GENE.SYMBOL','AveExpr','logFC','adj.P.Val','ENTREZ.ID','P.Value') -> diff.table
        
      top.t = diff.table[order(-diff.table$logFC),]
      st = top.t$logFC
      names(st) = top.t$GENE.SYMBOL
      fgseaRes <- fgsea(pathways = gset, nproc = 4, 
                        stats = st,
                        minSize=5,
                        maxSize=500,
                        nperm=10000)
     
      if(fgseaRes[1,'padj']<=0.1){
        
        print(
          plotEnrichment(gset[['P53.NEG.INT']],
                       st) + labs(title=paste(cancer,m,"NEG p53 int"))
        )
      }
      
      if(is.null(res)){
        res = fgseaRes[order(fgseaRes$padj),]
        res$contrast = m
        res$cancer = cancer
      }else{
        
        res.aux = fgseaRes[order(fgseaRes$padj),]
        res.aux$contrast = m
        res.aux$cancer = cancer
        res = rbind(res,res.aux)
        
      }
    }
  }
}
dev.off()
fwrite(res,paste(output.dir,'GSEA/gsea.pval.table.NEG.FC.csv',sep = ""))











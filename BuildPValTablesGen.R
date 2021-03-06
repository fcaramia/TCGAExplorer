rm(list = ls())
library(dplyr)
library(plyr)
library(Matrix)
library(gridExtra)
library(gtable)
library(grid)
#library(xlsx)
library(combinat)
library(data.table)
source("ExpressionPlotsFunctions.R")
datasets = c("ESCA","HNSC","LUSC","BLCA","LIHC","STAD","LGG","COAD","PAAD","READ","SKCM","LUAD")
input.dir = "/home/fcaramia/Documents/PhD/GenderAnalysis/TCGA/Analysis/TCGAExpressionExplorerOutput/"
output.dir = "/home/fcaramia/Documents/PhD/GenderAnalysis/TCGA/Analysis/TCGAExpressionExplorerOutput/"
signature.file =  "~/Documents/PhD/GenderAnalysis/TP53_interactions/Candidates.2018.12.20.full.csv"
sig.lists = fread(signature.file)
print.images = F
genes = sig.lists$GeneSymbol

###Start Main Loop###
print('Main Loop')
exp.mat = NULL
for (cancer in datasets)
{
  print(cancer)
  ###Create directories####
  #dir.create(paste(output.dir,cancer,sep = ""),showWarnings = F)
  dir.create(paste(output.dir,cancer,"/PValTables/",sep = ""),showWarnings = F)
  dir.create(paste(output.dir,cancer,"/CombinedExpressionTable/",sep = ""),showWarnings = F)
  #########################
  files = list.files(path = paste(input.dir,cancer,'/DifferentialExpression/',sep = ''), 
                     pattern = '.csv', all.files = FALSE,
             full.names = FALSE, recursive = FALSE,
             ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
  
  cancer.mat = NULL
  contrast.matrix = NULL
  pse = NULL
  for (f in files){

    if(file.exists(paste(input.dir,cancer,'/DifferentialExpression/',f,sep = ''))){
      
      diffexp = read.csv(paste(input.dir,cancer,'/DifferentialExpression/',f,sep = ''), as.is = T)
      diffexp %>% 
        select_('GENE.SYMBOL','adj.P.Val') %>%
        filter(GENE.SYMBOL%in%genes) -> lmatrix
      m = gsub("(.*)\\.csv",'\\1',f)
      names(lmatrix)[2] <- m
      pse = unique(c(pse,unlist(strsplit(m,'-'))))
      if (is.null(contrast.matrix)){
        contrast.matrix = lmatrix
      }else{
        lmatrix %>% 
          left_join(contrast.matrix, by = c("GENE.SYMBOL")) -> contrast.matrix
        
      }
      
      
      diffexp %>% 
        select_('GENE.SYMBOL','AveExpr','logFC','adj.P.Val') %>%
        filter(GENE.SYMBOL%in%genes) -> cmatrix
      
      if (is.null(cancer.mat)){
        colnames(cmatrix) = c('GENE.SYMBOL',
                                 lapply(X= colnames(cmatrix[,2:4]),
                                        FUN = function(x) return(paste(x,':',m,sep = ''))))
        cancer.mat = cmatrix
       
      }else{
        colnames(cmatrix) = c('GENE.SYMBOL',
                              lapply(X= colnames(cmatrix[,2:4]),
                                     FUN = function(x) return(paste(x,':',m,sep = ''))))
        cmatrix %>% 
          left_join(cancer.mat, by = c("GENE.SYMBOL")) -> cancer.mat
      }
      
    }
    
  }
  clns = colnames(cancer.mat)
  cancer.mat$DATA.SET = cancer
  cancer.mat = cancer.mat[,c('DATA.SET',clns)]
  
  write.csv(cancer.mat,paste(output.dir,cancer,'/CombinedExpressionTable/',cancer,'.csv',sep = ''),row.names = F)
  
  # if(file.exists(paste(output.dir,'CombinedExpressionTable.xlsx',sep = ''))){
  #   print('YES')
  #   write.xlsx(cancer.mat,paste(output.dir,'CombinedExpressionTable.xlsx',sep = ''),
  #              sheetName = cancer,row.names = F,append = T)
  # }else{
  #   print('No')
  #   write.xlsx(x = cancer.mat,file = paste(output.dir,'CombinedExpressionTable.xlsx',sep = ''),
  #              sheetName = cancer,row.names = F)
  # }
  if (print.images == T){
    for(g in contrast.matrix$GENE.SYMBOL){
      pic = as.data.frame(matrix(data=NA,nrow = length(pse),ncol = length(pse)))
      rownames(pic)= colnames(pic) = pse
      psen = length(pse)
      for (p1 in seq(1,psen-1)){
        for (p2 in seq(p1+1,psen)){
          
          ind = paste(pse[p1],pse[p2],sep = '-')
          if(ind %in% colnames(contrast.matrix)){
            val = contrast.matrix[which(contrast.matrix$GENE.SYMBOL==g),ind]
            #pic[pse[p1],pse[p2]] = ifelse(val<0.0001,format(val, scientific=T, digits = 4),format(val, scientific=F, digits = 2))
            pic[pse[p1],pse[p2]] = round(val,digits = 2)*100
            pic[pse[p2],pse[p1]] = pic[pse[p1],pse[p2]]
            
          }
          ind = paste(pse[p2],pse[p1],sep = '-')
          if(ind %in% colnames(contrast.matrix)){
            val = contrast.matrix[which(contrast.matrix$GENE.SYMBOL==g),ind]
            #pic[pse[p1],pse[p2]] = ifelse(val<0.0001,format(val, scientific=T, digits = 4),format(val, scientific=F, digits = 2))
            pic[pse[p1],pse[p2]] = round(val,digits = 2)*100
            pic[pse[p2],pse[p1]] = pic[pse[p1],pse[p2]]
            
          }
          
        }
      }
      #pic = format(round(pic,20),nsmall = 2)
      #pic = forceSymmetric(as.matrix(pic),'U')
      #ori.names2 = ori.names[ori.names%in%pse]
      #pic = data.frame(pic)
      #pic = pic[,ori.names2]
      #pic = pic[ori.names2,]
      #Replace with name
      #colnames(pic) = sapply(X =colnames(pic),FUN = function(x) return(h.table.names[[x]]))
      #rownames(pic) = sapply(X =rownames(pic),FUN = function(x) return(h.table.names[[x]]))
      tt1 = ttheme_minimal(rowhead = list(fg_params=list(fontface=c('bold'))))
      jpeg(file=paste(output.dir,cancer,'/PValTables/',g,'.jpeg',sep = ""),
           width = 4000, height = 2000, units = "px", res=170)
      sig_ind = Find_element_less_than(data.matrix(pic),25)
      gtab = tableGrob(pic)
      if (nrow(sig_ind)>0){
        for(i_gtab in seq(1:nrow(sig_ind))){
          ind_gtab = Find_cell(gtab,sig_ind[i_gtab,][1]+1,sig_ind[i_gtab,][2]+1)
          #print(ind_gtab)
          gtab$grobs[ind_gtab][[1]][['gp']] <- gpar(fontsize=15,fontface='bold')
        }
      }
      #grid.draw(gtab, theme = tt1)
      grid.draw(gtab)
      dev.off()
    }  
    
  }
  
  if (is.null(exp.mat)){
    
    exp.mat = cancer.mat
    
  }else{
    exp.mat = rbind.fill(exp.mat,cancer.mat)
  }
  
}

write.table(exp.mat,paste(output.dir,'DiffExpressionTable.txt',sep = ''),row.names = F ,sep = '\t', quote = FALSE)


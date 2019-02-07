rm(list = ls())
library(dplyr)
library(plyr)
library(Matrix)
library(gridExtra)
library(gtable)
library(grid)
library(xlsx)
library(combinat)
source("ExpressionPlotsFunctions.R")
demo = read.csv("~/Documents/PhD/GenderAnalysis/TCGA/Analysis/all.TCGA.curated.clinical.csv", as.is = T)
demo$PATIENT_ID = toupper(demo$PATIENT_ID)

datasets = c("ESCA","HNSC","LUSC","BLCA","LIHC","STAD","LGG","COAD","PAAD","READ","SKCM","LUAD")
#datasets = c("LGG",'ESCA')
#vars = list('Normal','Wt.TP53','Mt.TP53')
#vars = list('Normal','Tumor')
vars = list('Tumor')

pivots = list("female",'male')
input.dir = "/home/fcaramia/Documents/PhD/GenderAnalysis/TCGA/Analysis/TCGAExpressionExplorerOutput/"
output.dir = "/home/fcaramia/Documents/PhD/GenderAnalysis/TCGA/Analysis/TCGAExpressionExplorerOutput/"
signature.file =  "~/Documents/PhD/GenderAnalysis/TP53_interactions/TP53.list.3.txt"
sig.lists = read.delim(signature.file, sep = " ")

#t.order = c('F normal','M normal','F wt TP53','M wt TP53','F mt TP53','M mt TP53')
#t.order = c('F wt TP53','M wt TP53','F mt TP53','M mt TP53')
#t.order = c('F normal','M normal','F Tumour','M Tumour')
t.order = c('F Tumour','M Tumour')

#ori.names = c('Normal.female','Normal.male','Wt.TP53.female','Wt.TP53.male','Mt.TP53.female','Mt.TP53.male')
#ori.names = c('Wt.TP53.female','Wt.TP53.male','Mt.TP53.female','Mt.TP53.male')
ori.names = c('Normal.female','Normal.male','Tumor.female','Tumor.male')
ori.names = c('Tumor.female','Tumor.male')


h.table.names = hash(keys=ori.names, values = t.order)

genes = c()
for (s in rownames(sig.lists)){
  #Add signature
  gen = unlist(strsplit(as.character(sig.lists[s,'GENES']),',')) 
  genes = unique(c(genes,gen))
  
}
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
  cols = c()
  for (l in vars){
    for (p in pivots){
     f = paste(l,'.',p,sep = '')
     cols = c(cols,f)   
      
    }
  } 
  cancer.mat = NULL
  contrast.matrix = NULL
  pse = NULL
  ps = levels(as.factor(cols))
  psn = length(ps)
  for (p1 in seq(1,psn-1)){
    for (p2 in seq(p1+1,psn)){
      m = paste(ps[p1],ps[p2],sep = '-')
      #print (m)
      if(file.exists(paste(input.dir,cancer,'/DifferentialExpression/',m,'.csv',sep = ''))){
        pse = unique(c(pse,ps[p1],ps[p2]))
        diffexp = read.csv(paste(input.dir,cancer,'/DifferentialExpression/',m,'.csv',sep = ''), as.is = T)
        diffexp %>% 
          select_('GENE.SYMBOL','adj.P.Val') %>%
          filter(GENE.SYMBOL%in%genes) -> lmatrix
        names(lmatrix)[2] <- m
        
        if (is.null(contrast.matrix)){
          contrast.matrix = lmatrix
        }else{
          lmatrix %>% 
            left_join(contrast.matrix, by = c("GENE.SYMBOL")) -> contrast.matrix
          
        }
        
        
        diffexp %>% 
          select_('GENE.SYMBOL','logFC','adj.P.Val') %>%
          filter(GENE.SYMBOL%in%genes) -> cmatrix
        
        if (is.null(cancer.mat)){
          colnames(cmatrix) = c('GENE.SYMBOL',
                                   lapply(X= colnames(cmatrix[,2:3]),
                                          FUN = function(x) return(paste(x,':',m,sep = ''))))
          cancer.mat = cmatrix
         
        }else{
          colnames(cmatrix) = c('GENE.SYMBOL',
                                lapply(X= colnames(cmatrix[,2:3]),
                                       FUN = function(x) return(paste(x,':',m,sep = ''))))
          cmatrix %>% 
            left_join(cancer.mat, by = c("GENE.SYMBOL")) -> cancer.mat
        }
        
      }
    }
  }
  clns = colnames(cancer.mat)
  cancer.mat$DATA.SET = cancer
  cancer.mat = cancer.mat[,c('DATA.SET',clns)]
  
  ##Add numbers
  norm.counts = read.csv(paste(input.dir,cancer,'/Normalisation/Norm.Log.Counts.csv',sep = ""), as.is = T)
  rownames(norm.counts) = norm.counts[,1]
  norm.counts = norm.counts[,-1]
  
  PATIENT_IDs = gsub("TCGA\\.([[:alnum:]]{2})\\.([[:alnum:]]{4}).*","\\2", colnames(norm.counts))
  PATIENT_GENDER = demo[which(demo$PATIENT_ID%in%PATIENT_IDs),'GENDER']
  
  gt = table(PATIENT_GENDER)
  nf = gt['female']
  nm = gt['male']
  cancer.mat$'Females - Males' = paste(nf,nm,sep=' - ') 
  
  
  #write.csv(cancer.mat,paste(output.dir,cancer,'/CombinedExpressionTable/',cancer,'.csv',sep = ''),row.names = F)
  
  # if(file.exists(paste(output.dir,'CombinedExpressionTable.xlsx',sep = ''))){
  #   print('YES')
  #   write.xlsx(cancer.mat,paste(output.dir,'CombinedExpressionTable.xlsx',sep = ''),
  #              sheetName = cancer,row.names = F,append = T)
  # }else{
  #   print('No')
  #   write.xlsx(x = cancer.mat,file = paste(output.dir,'CombinedExpressionTable.xlsx',sep = ''),
  #              sheetName = cancer,row.names = F)
  # }
  for(g in contrast.matrix$GENE.SYMBOL){
    pic = as.data.frame(matrix(data=NA,nrow = length(pse),ncol = length(pse)))
    rownames(pic)= colnames(pic) = pse
    psen = length(pse)
    for (p1 in seq(1,psen-1)){
      for (p2 in seq(p1+1,psen)){
     
        ind = paste(pse[p1],pse[p2],sep = '-')
        if(ind %in% colnames(contrast.matrix))
          val = contrast.matrix[which(contrast.matrix$GENE.SYMBOL==g),ind]
          #pic[pse[p1],pse[p2]] = ifelse(val<0.0001,format(val, scientific=T, digits = 4),format(val, scientific=F, digits = 2))
          pic[pse[p1],pse[p2]] = round(val,digits = 2)*100
          pic[pse[p2],pse[p1]] = pic[pse[p1],pse[p2]]
          
      }
    }
    #pic = format(round(pic,20),nsmall = 2)
    #pic = forceSymmetric(as.matrix(pic),'U')
    ori.names2 = ori.names[ori.names%in%pse]
    pic = data.frame(pic)
    pic = pic[,ori.names2]
    pic = pic[ori.names2,]
    #Replace with name
    colnames(pic) = sapply(X =colnames(pic),FUN = function(x) return(h.table.names[[x]]))
    rownames(pic) = sapply(X =rownames(pic),FUN = function(x) return(h.table.names[[x]]))
    tt1 = ttheme_minimal(rowhead = list(fg_params=list(fontface=c('bold'))))
    jpeg(file=paste(output.dir,cancer,'/PValTables/',g,'_WtvMt.jpeg',sep = ""),
         width = 1500, height = 350, units = "px", res=170)
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
  if (is.null(exp.mat)){
    
    exp.mat = cancer.mat
    
  }else{
    exp.mat = rbind.fill(exp.mat,cancer.mat)
  }
  
}


#####special case
genomes.1000 = read.csv("/home/fcaramia/Documents/PhD/GenderAnalysis/1000.Genomes.Analysis/Output/DifferentialExpression/female-male.csv")
genomes.1000 = genomes.1000[which(genomes.1000$SYMBOL%in%genes),c('SYMBOL','logFC','adj.P.Val')]
genomes.1000$DATA.SET = '1000 Genomes'
genomes.1000$"Females:Males" = '333 - 327'
genomes.1000 = genomes.1000[,c('DATA.SET','SYMBOL','logFC','adj.P.Val','Females:Males')]

colnames(genomes.1000)  =  colnames(exp.mat)

res = rbind(exp.mat,genomes.1000)

write.table(res,paste(output.dir,'DiffExpressionTable1000GenomesAndTumors.tsv',sep = ''),row.names = F ,sep = '\t', quote = FALSE)


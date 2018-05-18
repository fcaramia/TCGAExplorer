library(limma)
library(dplyr)

GenerateRandomGeneSetBySize <- function(original.set, genes.to.choose, data.base){
  res = c()
  
  for(g in original.set){
    #print('Gene')
    #print(g)
    if(g%in%data.base$ENTREZ.ID){
      
      db = NULL
      size.scale = 0.1
      while (is.null(db))
      {
      
        len <- data.base[which(data.base$ENTREZ.ID == g),'GENE.LENGTH']
        top.len = len + len*size.scale
        bot.len = len - len*size.scale
        genes.to.choose %>% 
          filter(GENE.LENGTH>=bot.len & GENE.LENGTH <= top.len &
                   !(ENTREZ.ID%in%res) & !(ENTREZ.ID%in%original.set)) -> db
        size.scale = size.scale + 0.1
      }
      
      #print('Genes that match size:')
      #print(dim(db)[1])
      r.gene = sample(as.vector(db$ENTREZ.ID),1)
      res = c(r.gene,res)
    }else{
      print('Gene not in DB')
      print(g)
      genes.to.choose %>% 
        filter(!(ENTREZ.ID%in%res) & !(ENTREZ.ID%in%original.set)) -> db
      r.gene = sample(as.vector(db$ENTREZ.ID),1)
      res = c(r.gene,res)
    }
  }
  return(res)
}


DoDifferentialExpression <- function(mat,design,contrast,weights=NULL){
  
  #print('Fitting Linear Model with weights')
  fit <- lmFit(mat,design = design, weights = weights)
  #print('Fitting Contrast')
  fit2 = contrasts.fit(fit,contrast)
  fit2 = eBayes(fit2, trend = T, robust = T)
  rm(fit,mat,design,contrast,weights)
  return(fit2)
  
}

GetDiffExpTable <- function(contrast, fit, contrast.name, gene.annot){
  top.t = topTable(fit,coef = contrast.name, sort.by = 'p', number = Inf, p.value = 1,adjust.method = 'BH')
  top.t$ENTREZ.ID = as.integer(rownames(top.t))
  top.t$GENE.SYMBOL = gene.annot[match(rownames(top.t),gene.annot[,'Entrez.Gene.ID']),'Approved.Symbol']
  top.t$CHROMOSOME = gene.annot[match(rownames(top.t),gene.annot[,'Entrez.Gene.ID']),'Chrm']
  rm(contrast,fit,contrast.name,gene.annot)
  return(top.t)
}

GetAvgRanking <- function(mat, genes, col.rank, col.genes){
  s.mat = order(-abs(mat[,col.rank]))
  ind.genes = mat[,col.genes]%in%genes
  ret = mean(s.mat[ind.genes])
  rm(s.mat,ind.genes,mat,genes,col.rank,col.genes)
  return(ret)
  
}
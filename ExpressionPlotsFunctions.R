library(plyr)
library(reshape2)
library(ggplot2)
library(hash)
library(limma)
library(edgeR)
library(NanoStringNorm)
library(dplyr)
library(lazyeval)

Find_element_less_than <- function(table, value){
  inds = which(table <= value,arr.ind = T)
  return(inds)
}

Find_cell <- function(table,row,col, name='core-fg'){
  l <- table$layout
  which(l$t==row & l$l==col & l$name==name)
}

PrintDiffExpTable <- function(contrast, fit, dir, contrast.name, gene.annot){
  top.t = topTable(fit,coef = contrast.name, sort.by = 'p', number = Inf, p.value = 1,adjust.method = 'BH')
  top.t$ENTREZ.ID = as.integer(rownames(top.t))
  top.t$GENE.SYMBOL = gene.annot[match(rownames(top.t),gene.annot[,'Entrez.Gene.ID']),'Approved.Symbol']
  top.t$CHROMOSOME = gene.annot[match(rownames(top.t),gene.annot[,'Entrez.Gene.ID']),'Chrm']
  
  write.csv(top.t, paste(dir,contrast.name,'.csv',sep = ''), row.names = F)
  return(top.t)
}

PrintDiffExpTable2 <- function(contrast, fit, dir, contrast.name, gene.annot){
  top.t = topTable(fit,coef = contrast.name, sort.by = 'p', number = Inf, p.value = 1,adjust.method = 'BH')
  top.t$SYMBOL = rownames(top.t)
  top.t$SYMBOL.EXCEL.VIEW = paste('`',rownames(top.t),sep = '')
  top.t$CHROMOSOME = gene.annot[match(rownames(top.t),gene.annot[,'Approved.Symbol']),'Chrm']
  
  write.csv(top.t, paste(dir,contrast.name,'.csv',sep = ''), row.names = F)
  return(top.t)
}

DoDensityPlotGenderRMAF = function(dat, dat.col, fill.col, title.txt){
  
  print(
    ggplot(dat,aes_string(x = dat.col, colour= fill.col)) + geom_density() +
      scale_color_brewer(palette = 'Set1') + ggtitle(label = title.txt) + theme_minimal()
  ) 
}
DoDensityPlotFacetGenderRMAF = function(dat, dat.col,facet.col, fill.col, title.txt){
  fw = as.formula(paste("~", facet.col))
  print(
    ggplot(dat,aes_string(x = dat.col, colour= fill.col)) + geom_density() +
      scale_color_brewer(palette = 'Set1') + ggtitle(label = title.txt) +
      facet_wrap(fw, labeller = label_context) + theme_minimal()
  ) 
}

DoSmoothDualHistogram = function(dat, dat.col, fill.col, title.txt){
  
  select_(dat,dat.col,fill.col)%>%
    group_by_(fill.col)%>%
    dplyr::summarise_(med=interp(~median(var),var = as.name(dat.col))) -> cdat
  fg = as.formula(paste(fill.col,'~.',sep = ''))
  print(
    ggplot(dat,aes_string(x = dat.col, fill= fill.col)) +
      geom_histogram(bins = 15, alpha = 0.8, position = 'dodge') + 
      geom_vline(data=cdat, aes(xintercept=med,  colour=GENDER.TXT),
                                             linetype="dashed", size=1)+
      scale_fill_manual(values=c("pink1", "lightskyblue")) +ggtitle(label = title.txt) + theme_minimal()
  ) 
}

DoMutNumPlotGender = function(data,data.col,fill.col,x.lab,y.lab,title.txt, nmales, nfemales){
 
  data = as.data.frame(data)
 
  data[,fill.col] = as.factor(data[,fill.col])
  print(
      ggplot(data,aes_string(y=data.col,x=fill.col, fill = fill.col)) +
      stat_boxplot(geom ='errorbar',width = 0.2) + 
      geom_boxplot(outlier.colour = NA) + 
      geom_point(position = position_jitter(width = 0.2),show.legend = F, aes_string(shape = fill.col),alpha = 0.05) +
      labs(x=x.lab,y=y.lab) + ggtitle(label = title.txt)  + 
      scale_fill_manual(values=c("pink1", "lightskyblue")) + theme_minimal() +
      theme(legend.position = 'none',axis.ticks.x=element_blank(),
            axis.text.x = element_text(color='black'),strip.text = element_blank(),
            axis.title = element_text(color='black'), 
            plot.title = element_text(color='black',hjust = 0.5),
            axis.text.y = element_text(color='black')) +
        scale_x_discrete(labels = c('female'=paste('females',nfemales),'male'=paste('males',nmales)))
  )
  
}

DoMutNumPlotGenderFacet = function(data,data.col,fill.col,facet.col ,x.lab,y.lab,title.txt){
  
  
  fw = as.formula(paste("~", facet.col))
  print(
  ggplot(data,aes_string(y=data.col,x=fill.col, fill = fill.col)) +
    stat_boxplot(geom ='errorbar',width = 0.2) + 
    geom_boxplot(outlier.colour = NA) + 
    #geom_point(position = position_jitter(width = 0.2),show.legend = F, aes_string(shape = fill.col),alpha = 0.2) +
    labs(x=x.lab,y=y.lab) + ggtitle(label = title.txt) + facet_wrap(fw) +
    scale_fill_manual(values=c("pink1", "lightskyblue")) + theme_minimal() + coord_cartesian(ylim=c(0,5))
  )
  
}

DoMutNumViolinPlotGenderFacet = function(data,data.col,fill.col,facet.col ,x.lab,y.lab,title.txt, show.points = FALSE){
  
  
  fw = as.formula(paste("~", facet.col))
  plot = ggplot(data,aes_string(y=data.col,x=fill.col, fill = fill.col)) +
    stat_boxplot(geom ='errorbar',width = 0.2) + 
    geom_violin(outlier.colour = NA) + 
    labs(x=x.lab,y=y.lab) + ggtitle(label = title.txt) + facet_wrap(fw) +
    scale_fill_manual(values=c("pink1", "lightskyblue")) + theme_minimal() #+ coord_cartesian(ylim=c(0,1))
  
  if (show.points == TRUE){
    
    plot = plot + geom_point(position = position_jitter(width = 0.2),show.legend = F, 
                             aes_string(shape = fill.col),alpha = 0.1) 
    
  }
  
  print(
    plot
  )
  
}

DoMutNumViolinPlotGender = function(data,data.col,fill.col,facet.col ,x.lab,y.lab,title.txt, 
                                    show.points = FALSE, points.alpha = 0.01, point.col = NULL){
  
  plot = ggplot(data,aes_string(y=data.col,x=fill.col, fill = fill.col)) +
    stat_boxplot(geom ='errorbar',width = 0.2) + 
    geom_violin(outlier.colour = NA) + 
    labs(x=x.lab,y=y.lab) + ggtitle(label = title.txt)  +
    scale_fill_manual(values=c("pink1", "lightskyblue")) + theme_minimal() 
  
  if (show.points == TRUE){
    
    plot = plot + geom_point(position = position_jitter(width = 0.1), 
                             show.legend = F, aes_string(color = point.col),
                             alpha = points.alpha) + scale_color_manual(values = c('#000000FF','#FF0000FF')) 
      
  }
  
  print(
    plot
  )
  
}

  
DoMultBoxPlot = function(plot.data,cols.merge,cols.to.wrap,box.order,fill.color.col,ylab,xlab, box.title, 
                         plot.tittle, plot.type, var.colors = NULL)
{
  
  if (plot.type=='PDF')
    pdf(paste(plot.tittle,'.pdf',sep = ''))
  else
    jpeg(file=paste(plot.tittle,'.jpeg',sep = ""),width = 1200, height = 1400, units = "px", res=300)
    
  aux = plot.data
  if (length(cols.merge)>1){
    aux$MERGED <- apply( aux[ , cols.merge ] , 1 , paste , collapse = "." )
  }else{
    aux$MERGED <- aux[ , cols.merge ]
  }
  
  aux[,fill.color.col] = as.vector(aux[,fill.color.col])
  #for (i in unique(aux[,fill.color.col])){
  #  n = nrow(aux[which(aux[,fill.color.col]==i),])  
  #  n = paste('(',n,')',sep='')
  #  aux[which(aux[,fill.color.col]==i),fill.color.col] = paste(i,n)
  #}
  melted = melt(aux[,c('MERGED',cols.to.wrap,fill.color.col)])
  melted$MERGED = gsub(".*\\.(.*)","\\1",melted$MERGED)
  if (is.null(var.colors)){
    var.colors=rainbow(length(levels(as.factor(melted[,fill.color.col]))))  
  }
  
  ##Get numbers
  for (i in unique(melted$MERGED)){
    n = ''
    for (j in levels(as.factor(melted[,fill.color.col]))){
      
      n1 = nrow(melted[which(melted[,'MERGED']==i&melted[,fill.color.col]==j),])  
      #print(n1)
      n1=paste(n1,toupper(substr(j,1,1)))
      n = paste(n,n1,sep='')
    
    }
    n = paste('(',n,')',sep = '')
    melted[which(melted$MERGED==i),'MERGED'] = paste(i,n)
    if(!is.null(box.order)){
      box.order[which(box.order==i)] = paste(i,n)  
    }
    
    
  }
  if(!is.null(box.order))
  {
    melted$MERGED <- factor(melted$MERGED,levels = box.order)  
  }
  
   print(
     ggplot(melted,aes(x=MERGED,y=value, fill = melted[,fill.color.col])) + 
       stat_boxplot(geom ='errorbar') + geom_boxplot() +
       facet_wrap(~variable,labeller = ) + labs(y=ylab,x=xlab) + ggtitle(label = box.title) + 
       scale_fill_manual(fill.color.col,values=var.colors) + theme_minimal() + 
       theme(legend.position = 'none',axis.ticks.x=element_blank(),
             axis.text.x = element_text(color='black',face = 'bold'),strip.text = element_blank(),
             axis.title = element_text(color='black',face = 'bold'), 
             plot.title = element_text(color='black',face = 'bold',hjust = 0.5),
             axis.text.y = element_text(color='black',face = 'bold'))
   )
 dev.off()
   
}  
DoPValuesPlot <- function(plot.data,cols.merge,cols.to.wrap,ylab,xlab, box.title, plot.tittle){
  pdf(paste(plot.tittle,'.pdf',sep = ''))
  pd = plot.data
  n = length(cols.merge)
  for (i in seq(1,n-1)){
    for (j in seq(i+1,n)){
      col1 = cols.merge[i]
      col2 = cols.merge[j]
      for (v in cols.to.wrap){
        for(f in levels(as.factor(pd[,col1]))){
          ps = levels(as.factor(pd[,col2]))
          psn = length(ps)
          for (p1 in seq(1,psn-1)){
            for (p2 in seq(p1+1,psn)){
              pv1 = ps[p1]
              pv2 = ps[p2]
              d = pd[which(pd[,col1]==f&(pd[,col2]==pv1|pd[,col2]==pv2)),]
              pdd = melt(d[,c(col2,v)])
              d[,col2] = as.factor(d[,col2])
              p.val = 1
              pval = try(summary(lm(pdd[,'value']~pdd[,col2]))$coef[2,4])
              pv1.n = nrow(pd[which(pd[,col1]==f&pd[,col2]==pv1),])
              pv2.n = nrow(pd[which(pd[,col1]==f&pd[,col2]==pv2),])
              if (length(unique(d[,col2]))>1){
                #print(paste(pv1,pv2,f,sep = ' '))
                print(
                  ggplot(pdd,aes(x=col2,y=value, fill = pdd[,col2])) +
                    stat_boxplot(geom ='errorbar') + geom_boxplot() +
                    labs(y=ylab,x=col2) +  scale_fill_manual(col2, values = c('red','blue')) +
                    ggtitle(label = paste(v,', All',f,'by',col2,'\nLM p-val:',pval,'\n',pv1,':',pv1.n,';',pv2,":",pv2.n)) +
                    theme_minimal()
                )
              }
            }
          }
        }
        for(f in levels(as.factor(pd[,col2]))){
          ps = levels(as.factor(pd[,col1]))
          psn = length(ps)
          for (p1 in seq(1,psn-1)){
            for (p2 in seq(p1+1,psn)){
              pv1 = ps[p1]
              pv2 = ps[p2]
              d = pd[which(pd[,col2]==f&(pd[,col1]==pv1|pd[,col1]==pv2)),]
              pdd = melt(d[,c(col1,v)])
              d[,col1] = as.factor(d[,col1])
              p.val = 1
              pval = try(summary(lm(pdd[,'value']~pdd[,col1]))$coef[2,4])
              pv1.n = nrow(pd[which(pd[,col2]==f&pd[,col1]==pv1),])
              pv2.n = nrow(pd[which(pd[,col2]==f&pd[,col1]==pv2),])
              if (length(unique(d[,col1]))>1){
                #print(paste(pv1,pv2,f,sep = ' '))
                print(
                  ggplot(pdd,aes(x=col1,y=value, fill = pdd[,col1])) +
                    stat_boxplot(geom ='errorbar') + geom_boxplot() +
                    labs(y=ylab,x=col1) + scale_fill_manual(col1, values = c('red', 'blue')) +
                    ggtitle(label = paste(v,', All',f,'by',col1,'\nLM p-val:',pval,'\n',pv1,':',pv1.n,';',pv2,":",pv2.n)) +
                    theme_minimal()
                )
              }
            }
          }
        }
      }
    }
  }
  dev.off()  
}

DoNormPlots <-function(raw.data, norm.data, dir, tittle, annot, color.cols) {
  
  pdf(paste(dir,'/Normalisation.pdf',sep = ''))
  ####Before and After norms####
  par(mfrow= c(2,2))
  layout(matrix(c(1,2,3,4), 2, 2, byrow = TRUE))
  plot(density(log2(raw.data)),main = paste(tittle,'Raw',sep = '-'))
  plot(density(norm.data),main = paste(tittle,'Norm',sep = '-'))
  ####Boxplots#####
  boxplot(log2(raw.data[,1:50]), xlab="", ylab="Log2 Intensities",las=2,main="Unnormalised logIntensities")
  ## Let's add a blue horizontal line that corresponds to the median logCPM
  abline(h=median(log2(raw.data[,1:50])),col="blue")
  
  boxplot(norm.data[,1:50], xlab="", ylab="Log2 Intensities",las=2,main="Normalised logIntensities")
  ## Let's add a blue horizontal line that corresponds to the median logCPM
  abline(h=median(norm.data[,1:50]),col="blue")
  ###Norm Plots####
  meanSdPlot(data.matrix(norm.data))
  dev.off()
  
}

DoMDSPlot <- function(norm.data, dir, tittle, annot, color.cols){
  if(is.null(ncol(norm.data)))return()
  if(ncol(norm.data)<3)return()
  
  par(mfrow=c(1,1))
  for(i in color.cols){
    annot$color = 0
    col = 1
    for (j in levels(as.factor(annot[,i]))){
      annot[which(annot[,i]==j),'color'] = col
      col = col+1
    }
    ###MDS Plots###
    plotMDS(norm.data,main=paste(tittle,i,sep = '-'),
            labels =  annot[,i],
            col=annot[match(colnames(norm.data),annot$SAMPLE_ID),'color']
            ,gene.selection="pairwise",cex = 0.8)
    legend("topleft", legend=levels(as.factor(annot[,i])), col=levels(as.factor(annot$color)), pch=15)
    
  }
  
}

DoMultMDSPlots <- function(norm.data, dir, tittle, annot, color.cols){
  pdf(paste(dir,'/MDSPlots.pdf',sep = ''))
  print('Plot General MDS')
  DoMDSPlot(norm.data = norm.data,dir = dir, 
            tittle = tittle,annot = annot,color.cols = color.cols )
  n = length(color.cols)
  if(n==1){
    dev.off()
    return()
  }
  print('Plot detailed MDS')
  for (i in seq(1,n-1)){
    for (j in seq(i+1,n)){
      col1 = color.cols[i]
      col2 = color.cols[j]
      print(col1)
      print(col2)
      for(v in levels(as.factor(annot[,col1]))){
        samples = as.vector(annot[which(annot[,col1]==v),"SAMPLE_ID"])
        samp_annot = annot[which(annot[,col1]==v),]
        dat = norm.data[,samples]
        DoMDSPlot(norm.data = dat,dir = dir, 
                  tittle = paste(tittle,'All',v),annot = samp_annot,color.cols = c(col2) )
      }
      for(v in levels(as.factor(annot[,col2]))){
        samples =  as.vector(annot[which(annot[,col2]==v),"SAMPLE_ID"])
        dat = norm.data[,samples]
        samp_annot = annot[which(annot[,col2]==v),]
        DoMDSPlot(norm.data = dat,dir = dir, 
                  tittle = paste(tittle,'All',v),annot = samp_annot,color.cols = c(col1) )
      }
      
    }
  }
  dev.off()
  
}

DoExpressionPlots <- function(annot,i,comp.clinical.vars,comp.other.vars,comp.gene.muts,
                              h.gene.mut.order,sig.lists,h.other.var.order){
  
  for (c in comp.clinical.vars){
    annot.aux = annot
    annot = annot[which(!is.na(annot[,c])),]
    for (g in comp.gene.muts){
      plot.order = unlist(h.gene.mut.order[[g]])
      g2 = paste(g,'_STATUS',sep = '')
      annot = annot[which(!is.na(annot[,g2])),]
      for (s in rownames(sig.lists)){
        sig = as.character(sig.lists[s,'SIGNATURE'])
        #Plot signature alone
        plot.tittle = paste(output.dir,i,"/ExpressionPlots/",i,'.',sig,'.',g2,'.',c,sep = "") 
        
        DoMultBoxPlot(plot.data = annot, cols.merge = c(c,g2), cols.to.wrap = sig,
                      box.order = plot.order, fill.color.col = c,ylab = 'Log2 RNA counts',
                      xlab = g2,box.title = paste(i,sig,sep = '-'), plot.tittle, 'JPEG')
        # #Plot siganture with genes
        genes = unlist(strsplit(as.character(sig.lists[s,'GENES']),',')) 
        genes = genes[which(genes %in% colnames(annot))]
        # 
        # plot.tittle = paste(output.dir,i,"/ExpressionPlots/",i,'.',sig,'.with.genes.',g2,'.',c,sep = "") 
        # DoMultBoxPlot(plot.data = annot, cols.merge = c(c,g2), cols.to.wrap = c(sig,genes),
        #               box.order = plot.order, fill.color.col = c,ylab = 'Log2 RNA counts',
        #               xlab = g2,box.title = paste(i,sig,sep = '-'),plot.tittle, 'JPEG')
        
        # plot.tittle = paste(output.dir,i,"/ExpressionPlots/",i,'.',sig,'.all.pvals.',g2,'.',c,sep = "") 
        # DoPValuesPlot(plot.data = annot, cols.merge = c(c,g2), cols.to.wrap = c(sig,genes),
        #               ylab = 'RNA Log Fold Change',xlab = g2,
        #               box.title = paste(i,sig,sep = '-'), plot.tittle)
        
        #Plot genes individually
        for (gs in genes){
          if (gs%in%colnames(annot)){
            plot.tittle = paste(output.dir,i,"/ExpressionPlots/",i,'.',gs,'.',g2,'.',c,sep = "") 
            DoMultBoxPlot(plot.data = annot, cols.merge = c(c,g2), cols.to.wrap = gs,
                          box.order = plot.order, fill.color.col = c,ylab = 'Log2 RNA counts',
                          xlab = g2,box.title = paste(i,gs,sep = '-'), plot.tittle, 'JPEG')
            
          } 
        }
      }
      annot = annot.aux
    }
    
    for (v in comp.other.vars){
      plot.order = unlist(h.other.var.order[[v]])
      annot.aux = annot
      annot = annot[which(!is.na(annot[,c])),]
      for (s in rownames(sig.lists)){
        #Plot signature alone
        sig = as.character(sig.lists[s,'SIGNATURE'])
        plot.tittle = paste(output.dir,i,"/ExpressionPlots/",i,'.',sig,'.',v,'.',c,sep = "") 
        DoMultBoxPlot(plot.data = annot, cols.merge = c(c,v), cols.to.wrap = sig,
                      box.order = plot.order, fill.color.col = c,ylab = 'Log2 RNA counts',
                      xlab = paste(v),box.title = paste(i,sig,sep = '-'),plot.tittle, 'JPEG')
        
        #Plot siganture with genes
        genes = unlist(strsplit(as.character(sig.lists[s,'GENES']),',')) 
        genes = genes[which(genes %in% colnames(annot))]
        # plot.tittle = paste(output.dir,i,"/ExpressionPlots/",i,'.',sig,'.with.genes.',v,'.',c,sep = "") 
        # DoMultBoxPlot(plot.data = annot, cols.merge = c(c,v), cols.to.wrap = c(sig,genes),
        #               box.order = plot.order, fill.color.col = c,ylab = 'Log2 RNA counts',
        #               xlab = paste(v),box.title = paste(i,sig,sep = '-'), plot.tittle, 'JPEG')
         # plot.tittle = paste(output.dir,i,"/ExpressionPlots/",i,'.',sig,'.all.pvals.',v,'.',c,sep = "") 
         # DoPValuesPlot(plot.data = annot, cols.merge = c(c,v), cols.to.wrap = c(sig,genes),
         #               ylab = 'RNA Log Fold Change',xlab = paste(v),
         #               box.title = paste(i,sig,sep = '-'), plot.tittle)
         # 
        
        #Plot genes individually
        for (gs in genes){
          if (gs%in%colnames(annot)){
            plot.tittle = paste(output.dir,i,"/ExpressionPlots/",i,'.',gs,'.',v,'.',c,sep = "") 
            DoMultBoxPlot(plot.data = annot, cols.merge = c(c,v), cols.to.wrap = gs,
                          box.order = plot.order, fill.color.col = c,ylab = 'Log2 RNA counts',
                          xlab = paste(v),box.title = paste(i,gs,sep = '-'), plot.tittle, 'JPEG')
            
          }
        }
      }
      
    }
    annot = annot.aux
    
    rm(annot.aux,gs,sig,genes,v,s)
  }
}

DoCPMNorm <- function(raw.counts){
  
  ###Filter low read genes####
  L = min(colSums(raw.counts)) # min library size
  P = round(ncol(raw.counts)*0.50) #20% population
  dge = DGEList(counts = raw.counts)
  keep <- rowSums(cpm(dge) > 5/L*1e6) > P
  dge = DGEList(dge[keep,,keep.lib.sizes=F]) 
  ###########################
  ###Normalise####
  dge = calcNormFactors(dge)
  logCPM = cpm(dge,prior.count = 3, log = T)
  ###############
  
  return(logCPM)
  
}

NormaliseVoom <- function(raw.counts, design) {
  if (is.null(raw.counts)) {
    return()
  }
  ###Filter low read genes####
  L = min(colSums(raw.counts)) # min library size
  P = round(ncol(raw.counts)*0.20) #20% population
  dge = DGEList(counts = raw.counts)
  keep <- rowSums(cpm(dge) > 5/L*1e6) > P
  dge = DGEList(dge[keep,,keep.lib.sizes=F]) 
  ##Normalise
  dge = calcNormFactors(dge)
  v <- voom(dge, design, plot = T)
  return(v)
  
}

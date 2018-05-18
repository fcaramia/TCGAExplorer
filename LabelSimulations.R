rm(list = ls())
source("SimulationFunctions.R")
source("PropensityScoresFunctions.R")
library(dplyr)
no.simulations = 10
write.to.file = F
###Contrasts to test####
cont.compare = c('Mt.TP53.female-Mt.TP53.male','Wt.TP53.female-Wt.TP53.male','Tumor.female-Tumor.male')

###############
datasets = c("ESCA","HNSC","LUSC","BLCA","LIHC","STAD","LGG","COAD","PAAD","READ","SKCM")
datasets = c("LGG")
input.dir = "~/Documents/PhD/GenderAnalysis/TCGA/Analysis/TCGAExpressionExplorerOutput/"
output.dir = "~/Documents/PhD/GenderAnalysis/TCGA/Analysis/TCGAExpressionExplorerOutput/"
comp.gene.muts = c("TP53")
comp.other.vars = c("SAMPLE_TYPE")
comp.clinical.vars = c('GENDER') #PIVOT

####Read Annotation file####
gene.annot = read.csv("~/Documents/PhD/GenderAnalysis/TCGA/Analysis/annotation.all.genes.csv",as.is = T)
gene.annot = gene.annot[,-1]
gene.annot = gene.annot[gene.annot$Chrm!='',]
gene.annot = gene.annot[order(gene.annot$Chrm,gene.annot$Entrez.Gene.ID),]
#########################

#Cofounding Factors definitions 
cofounding.factors.ind = c('SMOKING_STATUS','RACE','PATHOLOGIC_STAGE')
cofounding.factors.all = c('SMOKING_STATUS','RACE')
cofounding.factors.default = c('PATIENT_ID','GENDER','AGE_AT_INITIAL_PATHOLOGIC_DIAGNOSIS')
##############################

###Mutations####
mutations = read.csv("~/Documents/PhD/GenderAnalysis/TCGA/Analysis/all.reduced.TCGA.curated.mutations.csv", as.is = T)
mutations$PATIENT_ID = gsub("TCGA\\-([[:alnum:]]{2})\\-([[:alnum:]]{4}).*","\\2", mutations$SAMPLE_ID)
###############

###Clinical Data####
demo = read.csv("~/Documents/PhD/GenderAnalysis/TCGA/Analysis/all.TCGA.curated.clinical.csv", as.is = T)
demo$PATIENT_ID = toupper(demo$PATIENT_ID)


signature.file =  "~/Documents/PhD/GenderAnalysis/TCGA/Analysis/TP53.lists.txt"
###Read Signatures####
sig.lists = read.delim(signature.file, sep = " ")

##Big loop
for(cancer in datasets){
  
  print(cancer)
  ###Create directories####
  dir.create(paste(output.dir,cancer,sep = ""),showWarnings = F)
  #dir.create(paste(output.dir,cancer,"/DifferentialExpression",sep = ""),showWarnings = F)
  #########################
  ####Clinical data for cancer
  clinical_data = filter(demo,CANCER_TYPE==cancer)
  cf = cofounding.factors.ind
  cfd = cofounding.factors.default
  ###Read Norm counts#####
  print('Reading Norm counts...')
  norm.counts = read.csv(paste(input.dir,cancer,'/Normalisation/Norm.Log.Counts.csv',sep = ""), as.is = T)
  rownames(norm.counts) = norm.counts[,1]
  norm.counts = norm.counts[,-1]
  #########################
  
  ###Only genes of interest####
  #norm.counts = norm.counts[which(rownames(norm.counts)%in%use.genes$Entrez.Gene.ID),]
  
  ####Prepare annotation#####
  annot = as.data.frame(cbind(SAMPLE_EXP=colnames(norm.counts)))
  annot$PATIENT_ID = gsub("TCGA\\.([[:alnum:]]{2})\\.([[:alnum:]]{4}).*","\\2", annot[,1])
  annot$SAMPLE_ID = gsub("TCGA\\.([[:alnum:]]{2})\\.([[:alnum:]]{4})\\.([[:alnum:]]{3}).*",
                         "TCGA\\-\\1\\-\\2\\-\\3", annot[,1])
  
  ##Use counts with clinical data
  keep = annot$PATIENT_ID%in%clinical_data$PATIENT_ID
  annot = annot[keep,]
  norm.counts = norm.counts[,annot$SAMPLE_EXP]
  #Add Clinical Variables and filter samples in weights
  annot %>%
    #left_join(demo[,c('PATIENT_ID',comp.clinical.vars)]) -> annot
    left_join(clinical_data[,unique(c('PATIENT_ID',comp.clinical.vars,cf,cfd))]) -> annot
  
  #Add Sample_Type#
  annot$SAMPLE_CODE = gsub("TCGA\\.([[:alnum:]]{2})\\.([[:alnum:]]{4}\\.([[:alnum:]]{2})).*","\\3", annot[,1])
  annot$SAMPLE_TYPE <- ifelse(as.integer(annot$SAMPLE_CODE)>=20, 'Control',
                              ifelse(as.integer(annot$SAMPLE_CODE)>=10, 'Normal', 'Tumor'))
  
  #IF multiple samples from the same patient, Use lowest code 
  annot %>% 
    group_by(PATIENT_ID,SAMPLE_TYPE) %>% dplyr::summarise(MIN.CODE = min(SAMPLE_CODE)) -> aux
  
  aux$ID = paste(aux$PATIENT_ID,aux$MIN.CODE,sep = '')
  annot$auxID = paste(annot$PATIENT_ID,annot$SAMPLE_CODE,sep = '')
  annot = filter(annot,auxID%in%aux$ID)
  annot$auxID = NULL
  rm(aux)
  
  
  #Add Mutants
  mutations %>% 
    filter(CANCER_TYPE == cancer) -> muts.dataset
  for (g in comp.gene.muts){
    
    g2 = paste(g,'_STATUS',sep = '')
    muts.dataset %>% 
      filter(HUGO_SYMBOL == g) -> muts.gene
    annot[,g2] = NA
    
    annot[which(annot$SAMPLE_TYPE == 'Normal'), g2] = 'Normal'
    annot[which(annot$PATIENT_ID %in% muts.gene$PATIENT_ID&annot$SAMPLE_TYPE=='Tumor'),g2] = paste('Mt.',g,sep = '')
    annot[which(annot$SAMPLE_TYPE=='Tumor'&is.na(annot[,g2])),g2] = paste('Wt.',g,sep = '')
  }
  ################################
  rm(muts.gene,g2)
  ####Build contrast matrix
  
  new_columns = list()
  for (g in comp.gene.muts){
    g2 = paste(g,'_STATUS',sep = '')
    
    for (cv in comp.clinical.vars){
      nc = paste(g2,cv,sep = '.')
      new_columns = c(new_columns,nc)
      annot[,nc] = paste(annot[,g2],annot[,cv],sep='.')
      
    }
    nc = paste(g2,'all',sep = '.')
    new_columns = c(new_columns,nc)
    annot[,nc] = annot[,g2]
  }
  
  for (g in comp.other.vars){
    for (cv in comp.clinical.vars){
      nc = paste(g,cv,sep = '.')
      annot[,nc] = paste(annot[,g],annot[,cv],sep='.')
      new_columns = c(new_columns,nc)
    }
    nc = paste(g,'all',sep = '.')
    new_columns = c(new_columns,nc)
    annot[,nc] = annot[,g]
  }
  #print(new_columns)
  #print("Building Contrasts and doing DE Analysis")
  for (ncs in new_columns)
  {
    vals = levels(factor(annot[,ncs]))
    n = length(vals)
    if(n>1){   
      design = model.matrix(~0+annot[,ncs])
      colnames(design) = vals
      rownames(design) = as.vector(annot$SAMPLE_EXP)
      counts = data.matrix(norm.counts[,as.vector(annot$SAMPLE_EXP)])
      for (i in seq(1,n-1)){
        for (j in seq(i+1,n)){
          col1 = vals[i]
          col2 = vals[j]
          cont = paste(col1,col2,sep = '-')
          
          if(cont %in% cont.compare){
            term1.1 = gsub('(.*)\\..*','\\1',col1)
            term2.1 = gsub('(.*)\\..*','\\1',col2)
            term1.2 = gsub('.*\\.(.*)','\\1',col1)
            term2.2 = gsub('.*\\.(.*)','\\1',col2)
            print(cont)
            des = design[,c(col1,col2)]
            keep = rowSums(des) > 0
            des = des[keep,]
            counts.small = data.matrix(counts[,rownames(des)])
            annot.small = filter(annot,SAMPLE_EXP %in% rownames(des))
            
            #Calculating cofounding factors
            print('Calculating cofounding factor weights')
            weights.small = NULL
            if(term1.2!=term2.2&&term1.1==term2.1&&term1.1!='Normal'){
              weights.small = DoCalcWeights(clinical.data = annot.small,default.confounding.factors = cfd,
                                            extra.cofounding.factors = cf)
              
              #Filter samples with weights
              if(!is.null(weights.small))
              {
                annot.small = filter(annot.small, PATIENT_ID%in% weights.small$PATIENT_ID)
                des = des[as.vector(annot.small$SAMPLE_EXP),]
                counts.small = counts.small[,rownames(des)]
              }
              
            }
            rm(keep)
            #print('Building Contrast')
            contrast = makeContrasts(paste(col1,col2,sep = '-'),levels = des)
            df = DoDifferentialExpression(mat=counts.small,design = des,contrast = contrast,
                                          weights = weights.small$WEIGHTS)
            #print('Printing Results')
            df.table = GetDiffExpTable(contrast = contrast, fit = df, 
                                       contrast.name = paste(col1,col2,sep = '-'), gene.annot)
            print(head(df.table))
            ori.score = data.frame(matrix(ncol=nrow(sig.lists),nrow = 1,data = NA))
            colnames(ori.score) = sig.lists$SIGNATURE
            for (s in rownames(sig.lists)){
              genes = unlist(strsplit(as.character(sig.lists[s,'GENES']),',')) 
              signature = as.vector(sig.lists[s,'SIGNATURE'])
              ori.score[1,signature]=GetAvgRanking(mat = df.table, genes = genes, col.rank = 't', col.genes = 'GENE.SYMBOL')
            }
            rm(df,df.table)
            
            ##Generate Simulations
            sim.list = vector(mode='character',length=no.simulations)
            scores = data.frame(matrix(ncol=nrow(sig.lists),nrow = no.simulations ,data = 0) )
            colnames(scores) = sig.lists$SIGNATURE
            for(s in seq(1:no.simulations)){
              prct = (s/no.simulations)*100
              if(prct%%10==0){
                print(paste(prct,'%',sep=''))
              }
              #set.seed(s)
              sim = sample(seq(0,1),size = nrow(des), replace = T)
              sim.str = paste(sim,sep = '',collapse = '')
              while (sim.str%in%sim.list){
                print('repeat')
                sim = sample(seq(0,1),size = nrow(des), replace = T)
                sim.str = paste(sim,sep = '',collapse = '')
              }
              sim.list[s] = sim.str
              des[,1] = sim
              des[,2] = 1-sim
              df.sim = DoDifferentialExpression(mat=counts.small,design = des,contrast = contrast,
                                                weights = NULL)
              
              df.table.sim = GetDiffExpTable(contrast = contrast, fit = df.sim, 
                                             contrast.name = paste(col1,col2,sep = '-'), gene.annot)
              #Check rank of signature
              for (sl in rownames(sig.lists)){
                genes = unlist(strsplit(as.character(sig.lists[sl,'GENES']),',')) 
                signature = as.vector(sig.lists[sl,'SIGNATURE'])
                #Perform task for simulation
                sim.score=GetAvgRanking(mat = df.table.sim, genes = genes, col.rank = 't', col.genes = 'GENE.SYMBOL')
                scores[s,signature]=sim.score
              }
              #Get rid of trash
              rm(df.sim,df.table.sim,sim,sim.str,sim.score)
              gc()
            }
            gc()
            ##Check simulation results
            for (sl in rownames(sig.lists)){
              signature = as.vector(sig.lists[sl,'SIGNATURE'])
              sim.pval = rank(as.numeric(append(c(ori.score[1,signature]),scores[,signature])))[1] / no.simulations
              #print(ori.score)
              #print(scores)
              print(paste(cancer,cont,signature,sim.pval,sep = ','))
              if(write.to.file==T){
                write(paste(cancer,cont,signature,sim.pval,sep = ','),
                      file = paste(output.dir,cont,'.',no.simulations,'.NewSimulations.csv',sep = '')
                      ,append = T)    
              }
              
            }
            
          }
          
        }
      }
      rm(design,counts,vals)
    }
  }
  
  ########################
  rm(norm.counts,annot)
}


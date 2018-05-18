rm(list = ls())
library(gtools)
source("SimulationFunctions.R")
source("ExpressionPlotsFunctions.R")
library(dplyr)
no.simulations = 10000
write.to.file = T
###Contrasts to test####
cont.compare = c('female,male')
###############
datasets = c("ESCA","HNSC","LUSC","BLCA","LIHC","STAD","LGG","COAD","PAAD","READ","SKCM")
input.dir = "~/Documents/PhD/GenderAnalysis/TCGA/Analysis/TCGAExpressionExplorerOutput/"
output.dir = "~/Documents/PhD/GenderAnalysis/TCGA/Analysis/TCGAExpressionExplorerOutput/"

##Read TP53 interactors###
tp53.auto = read.delim("~/Documents/PhD/GenderAnalysis/TP53_interactions/TP53.AUTO.Interactors.txt", sep = '\t')

########

###Gene lenghts

gene.lenghts = read.csv("~/Documents/PhD/GenderAnalysis/TCGA/Analysis/All.genes.lengths.csv")
all.gene.lenghts = gene.lenghts
gene.lenghts %>% filter(CHROMOSOME!='X' & CHROMOSOME!='Y' ) -> gene.lenghts.auto
gene.lenghts %>% filter(CHROMOSOME=='X') -> gene.lenghts.x
gene.lenghts %>% filter(GENE.SYMBOL%in%tp53.auto$GENES) -> gene.lenghts.tp53.auto


random.genes = gene.lenghts.auto

###Clinical Data####
demo = read.csv("~/Documents/PhD/GenderAnalysis/TCGA/Analysis/all.TCGA.curated.clinical.csv", as.is = T)
demo$PATIENT_ID = toupper(demo$PATIENT_ID)
##RMAF DATA
RMAF = read.csv("~/Documents/PhD/GenderAnalysis/TCGA/Analysis/all.TCGA.ExMut.csv", as.is=T)
RMAF %>% 
  left_join(demo[c('GENDER','PATIENT_ID')]) %>%
  filter(!is.na(GENDER)) %>% 
  group_by(PATIENT_ID) %>% 
  dplyr::summarise(N.MUTS=n()) %>% filter(N.MUTS<10000) -> validated.patients

RMAF %>% 
  left_join(demo[c('GENDER','PATIENT_ID')]) %>%
  filter(!is.na(GENDER)) %>% 
  filter(PATIENT_ID%in%validated.patients$PATIENT_ID) -> RMAF 




#Clean low reads
RMAF$SUM = RMAF$REF_COUNT + RMAF$ALT_COUNT
RMAF %>% filter(SUM>10) -> RMAF
RMAF$RMAF = RMAF$ALT_COUNT/RMAF$SUM
RMAF$LOGIT.RMAF = logit(x = RMAF$RMAF,min = 0-0.1 , max = 1+0.1)


##Do something about noise?######


signature.file =  "~/Documents/PhD/GenderAnalysis/TP53_interactions/TP53.AUTO.Signatures.txt"
###Read Signatures####
sig.lists = read.delim(signature.file, sep = " ")

RMAF$GENDER.TXT = RMAF$GENDER
RMAF$GENDER = ifelse(RMAF$GENDER=='male',1,0)

RMAF.bck = RMAF
# for (cancer in datasets)
# {
#   print(cancer)  
#   RMAF = RMAF.bck
#   RMAF %>% filter(CANCER_TYPE==cancer) -> RMAF
#   if(nrow(RMAF)<100)
#   {
#     print('Skip')
#     next
#   }
  #Compute original p-val
ori.score = data.frame(matrix(ncol=nrow(sig.lists),nrow = 1,data = NA))
colnames(ori.score) = sig.lists$SIGNATURE
for (s in rownames(sig.lists)){
  
  genes = unlist(strsplit(as.character(sig.lists[s,'GENES']),',')) 
  RMAF.S = filter(RMAF,HUGO_SYMBOL%in%genes)
  genes = unique(RMAF.S$HUGO_SYMBOL)
  signature = as.vector(sig.lists[s,'SIGNATURE'])
  #print(signature)
  pvalue <- try(summary(lm(data = RMAF.S,formula = LOGIT.RMAF~GENDER))$coef[2,4])
  pvalue.all = format(pvalue, nsmall = 2)
  print(pvalue.all)
  ori.score[1,signature] =pvalue.all  
}
#Compute simulations
sim.list = matrix(ncol=nrow(sig.lists),nrow = no.simulations ,data = "") 
scores = data.frame(matrix(ncol=nrow(sig.lists),nrow = no.simulations ,data = 0) )
colnames(scores) = sig.lists$SIGNATURE
colnames(sim.list) = sig.lists$SIGNATURE
for(s in seq(1:no.simulations)){
  prct = (s/no.simulations)*100
  print (s)
  if(prct%%10==0){
    print(paste(prct,'%',sep=''))
  }
  #Check rank of signature
  for (sl in rownames(sig.lists)){
    genes = unlist(strsplit(as.character(sig.lists[sl,'GENES']),',')) 
    RMAF.S = filter(RMAF,HUGO_SYMBOL%in%genes)
    genes = unique(RMAF.S$HUGO_SYMBOL)
    signature = as.vector(sig.lists[sl,'SIGNATURE'])
    n.genes = length(genes)
    genes.eID = all.gene.lenghts[which(genes%in%all.gene.lenghts$GENE.SYMBOL),"ENTREZ.ID"]
    #print(n.genes)
    db = random.genes[which(random.genes$GENE.SYMBOL%in%RMAF$HUGO_SYMBOL),]
    random.list = GenerateRandomGeneSetBySize(original.set = genes.eID,data.base = all.gene.lenghts, genes.to.choose = db)
    random.vector = sort(as.vector(random.list))
    sim.str = paste(random.vector,sep = '',collapse = '')
    #print(random.vector)
    while (sim.str%in%sim.list[,signature]){
      print('repeated sim')
      db = random.genes[which(random.genes$GENE.SYMBOL%in%RMAF$HUGO_SYMBOL),]
      random.list = GenerateRandomGeneSetBySize(original.set = genes.eID,data.base = all.gene.lenghts, genes.to.choose = db)
      random.vector = sort(as.vector(random.list))
      sim.str = paste(random.vector,sep = '',collapse = '')
    }
    sim.list[s,signature] = sim.str
    #print(length(random.vector))
    
    rg = random.genes[which(random.genes$ENTREZ.ID%in%random.vector),'GENE.SYMBOL']
    #Perform task for simulation
    sim.dat = as.data.frame(filter(RMAF,HUGO_SYMBOL%in%rg))
    pvalue = Inf
    pvalue <- try(summary(lm(data = sim.dat,formula = LOGIT.RMAF~GENDER))$coef[2,4])
    pvalue.all = format(pvalue, nsmall = 2)
    scores[s,signature]=pvalue.all
  }
}
##Check simulation results
for (sl in rownames(sig.lists)){
  signature = as.vector(sig.lists[sl,'SIGNATURE'])
  sim.pval = rank(as.numeric(append(c(ori.score[1,signature]),scores[,signature])))[1] / no.simulations
  #print(ori.score)
  #print(scores)
  print(paste('RMAF',signature,sim.pval,sep = ','))
  if(write.to.file==T){
    write(paste('ALL.datasets.AUTO.V2',signature,sim.pval,sep = ','),
          file = paste(output.dir,'RMAF','.',no.simulations,'.Simulations.csv',sep = '')
          ,append = T)    
  }
}

#}



# ##Big loop
# for(cancer in datasets){
#   
#   print(cancer)
#   ###Create directories####
#   dir.create(paste(output.dir,cancer,sep = ""),showWarnings = F)
#   dir.create(paste(output.dir,cancer,"/DifferentialExpression",sep = ""),showWarnings = F)
#   #########################
#   
#   ###Read Norm counts#####
#   print('Reading Norm counts...')
#   norm.counts = read.csv(paste(input.dir,cancer,'/Normalisation/Norm.Log.Counts.csv',sep = ""), as.is = T)
#   rownames(norm.counts) = norm.counts[,1]
#   norm.counts = norm.counts[,-1]
#   #########################
#   
#   ###Only genes of interest####
#   #norm.counts = norm.counts[which(rownames(norm.counts)%in%use.genes$Entrez.Gene.ID),]
#   
#   ####Prepare annotation#####
#   annot = as.data.frame(cbind(SAMPLE_EXP=colnames(norm.counts)))
#   annot$PATIENT_ID = gsub("TCGA\\.([[:alnum:]]{2})\\.([[:alnum:]]{4}).*","\\2", annot[,1])
#   annot$SAMPLE_ID = gsub("TCGA\\.([[:alnum:]]{2})\\.([[:alnum:]]{4})\\.([[:alnum:]]{3}).*",
#                          "TCGA\\-\\1\\-\\2\\-\\3", annot[,1])
#   
#   ###Propensity Scores per Tumor Type###
#   weights.tumor = read.csv(paste(output.dir,cancer,'/PS/PS.csv',sep = ""))
#   
#   #Add Clinical Variables and filter samples in weights
#   annot %>%
#     #left_join(demo[,c('PATIENT_ID',comp.clinical.vars)]) -> annot
#     left_join(demo[,c('PATIENT_ID',comp.clinical.vars)]) %>%
#     filter(PATIENT_ID %in% weights.tumor$PATIENT_ID) %>%
#     left_join(weights.tumor) -> annot
#   
#   #Add Sample_Type#
#   annot$SAMPLE_CODE = gsub("TCGA\\.([[:alnum:]]{2})\\.([[:alnum:]]{4}\\.([[:alnum:]]{2})).*","\\3", annot[,1])
#   annot$SAMPLE_TYPE <- ifelse(as.integer(annot$SAMPLE_CODE)>=20, 'Control',
#                               ifelse(as.integer(annot$SAMPLE_CODE)>=10, 'Normal', 'Tumor'))
#   
#   #Add Mutants
#   mutations %>% 
#     filter(CANCER_TYPE == cancer) -> muts.dataset
#   for (g in comp.gene.muts){
#     
#     g2 = paste(g,'_STATUS',sep = '')
#     muts.dataset %>% 
#       filter(HUGO_SYMBOL == g) -> muts.gene
#     annot[,g2] = NA
#     
#     annot[which(annot$SAMPLE_TYPE == 'Normal'), g2] = 'Normal'
#     annot[which(annot$PATIENT_ID %in% muts.gene$PATIENT_ID&annot$SAMPLE_TYPE=='Tumor'),g2] = paste('Mt.',g,sep = '')
#     annot[which(annot$SAMPLE_TYPE=='Tumor'&is.na(annot[,g2])),g2] = paste('Wt.',g,sep = '')
#   }
#   ################################
#   rm(muts.gene,g2)
#   ####Build contrast matrix
#   
#   new_columns = list()
#   for (g in comp.gene.muts){
#     g2 = paste(g,'_STATUS',sep = '')
#     
#     for (cv in comp.clinical.vars){
#       nc = paste(g2,cv,sep = '.')
#       new_columns = c(new_columns,nc)
#       annot[,nc] = paste(annot[,g2],annot[,cv],sep='.')
#       
#     }
#     nc = paste(g2,'all',sep = '.')
#     new_columns = c(new_columns,nc)
#     annot[,nc] = annot[,g2]
#   }
#   
#   for (g in comp.other.vars){
#     for (cv in comp.clinical.vars){
#       nc = paste(g,cv,sep = '.')
#       annot[,nc] = paste(annot[,g],annot[,cv],sep='.')
#       new_columns = c(new_columns,nc)
#     }
#     nc = paste(g,'all',sep = '.')
#     new_columns = c(new_columns,nc)
#     annot[,nc] = annot[,g]
#   }
#   #print(new_columns)
#   #print("Building Contrasts and doing DE Analysis")
#   for (ncs in new_columns)
#   {
#     vals = levels(factor(annot[,ncs]))
#     n = length(vals)
#     if(n>1){   
#       design = model.matrix(~0+annot[,ncs])
#       colnames(design) = vals
#       rownames(design) = annot$SAMPLE_EXP
#       counts = data.matrix(norm.counts[,annot$SAMPLE_EXP])
#       for (i in seq(1,n-1)){
#         for (j in seq(i+1,n)){
#           col1 = vals[i]
#           col2 = vals[j]
#           cont = paste(col1,col2,sep = '-')
#           
#           if(cont %in% cont.compare){
#             print(cont)
#             des = design[,c(col1,col2)]
#             keep = rowSums(des) > 0
#             des = des[keep,]
#             counts.small = counts[,keep]
#             annot.small = filter(annot,SAMPLE_EXP %in% rownames(des))
#             rm(keep)
#             #print('Building Contrast')
#             contrast = makeContrasts(paste(col1,col2,sep = '-'),levels = des)
#             df = DoDifferentialExpression(mat=counts.small,design = des,contrast = contrast,
#                                           weights = annot.small$WEIGHTS)
#             #print('Printing Results')
#             df.table = GetDiffExpTable(contrast = contrast, fit = df, 
#                                        contrast.name = paste(col1,col2,sep = '-'), gene.annot)
#             
#             ori.score = data.frame(matrix(ncol=nrow(sig.lists),nrow = 1,data = NA))
#             colnames(ori.score) = sig.lists$SIGNATURE
#             for (s in rownames(sig.lists)){
#               genes = unlist(strsplit(as.character(sig.lists[s,'GENES']),',')) 
#               signature = as.vector(sig.lists[s,'SIGNATURE'])
#               ori.score[1,signature]=GetAvgRanking(mat = df.table, genes = genes, col.rank = 't', col.genes = 'GENE.SYMBOL')
#             }
#             rm(df,df.table)
#             
#             ##Generate Simulations
#             sim.list = vector(mode='character',length=no.simulations)
#             scores = data.frame(matrix(ncol=nrow(sig.lists),nrow = no.simulations ,data = 0) )
#             colnames(scores) = sig.lists$SIGNATURE
#             for(s in seq(1:no.simulations)){
#               prct = (s/no.simulations)*100
#               if(prct%%10==0){
#                 print(paste(prct,'%',sep=''))
#               }
#               #set.seed(s)
#               sim = sample(seq(0,1),size = nrow(des), replace = T)
#               sim.str = paste(sim,sep = '',collapse = '')
#               while (sim.str%in%sim.list){
#                 print('repeat')
#                 sim = sample(seq(0,1),size = nrow(des), replace = T)
#                 sim.str = paste(sim,sep = '',collapse = '')
#               }
#               sim.list[s] = sim.str
#               des[,1] = sim
#               des[,2] = 1-sim
#               df.sim = DoDifferentialExpression(mat=counts.small,design = des,contrast = contrast,
#                                                 weights = annot.small$WEIGHTS)
#               
#               df.table.sim = GetDiffExpTable(contrast = contrast, fit = df.sim, 
#                                              contrast.name = paste(col1,col2,sep = '-'), gene.annot)
#               #Check rank of signature
#               for (sl in rownames(sig.lists)){
#                 genes = unlist(strsplit(as.character(sig.lists[sl,'GENES']),',')) 
#                 signature = as.vector(sig.lists[sl,'SIGNATURE'])
#                 #Perform task for simulation
#                 sim.score=GetAvgRanking(mat = df.table.sim, genes = genes, col.rank = 't', col.genes = 'GENE.SYMBOL')
#                 scores[s,signature]=sim.score
#               }
#               #Get rid of trash
#               rm(df.sim,df.table.sim,sim,sim.str,sim.score)
#               gc()
#             }
#             gc()
#             ##Check simulation results
#             for (sl in rownames(sig.lists)){
#               signature = as.vector(sig.lists[sl,'SIGNATURE'])
#               sim.pval = rank(append(c(ori.score[1,signature]),scores[,signature]))[1] / no.simulations
#               #print(ori.score)
#               #print(scores)
#               print(paste(cancer,cont,signature,sim.pval,sep = ','))
#               if(write.to.file==T){
#                 write(paste(cancer,cont,signature,sim.pval,sep = ','),
#                       file = paste(output.dir,cont,'.',no.simulations,'.Simulations.csv',sep = '')
#                       ,append = T)    
#               }
#               
#             }
#             
#           }
#           
#         }
#       }
#       rm(design,counts,vals)
#     }
#   }
#   
#   ########################
#   rm(norm.counts,annot,weights.tumor)
# }
# 
# 
# RMAF %>% group_by(HUGO_SYMBOL,GENDER.TXT) %>% dplyr::summarise(MED.RMAF = median(LOGIT.RMAF),N = n()) -> l
# res = matrix(nrow = length(unique(l$HUGO_SYMBOL)),ncol = 3 , data = 0)
# rownames(res) = unique(l$HUGO_SYMBOL)
# for(g in unique(l$HUGO_SYMBOL))
# {
#   filter(l,HUGO_SYMBOL==g&GENDER.TXT=='female') -> f
#   filter(l,HUGO_SYMBOL==g&GENDER.TXT=='male') -> m
# 
#   if (dim(m)[1]!=0 & dim(f)[1]!=0)
#   {
# 
#     res[g,1] = g
#     res[g,2] = abs(f$MED.RMAF-m$MED.RMAF)
#     res[g,3] = min(f$N,m$N)
#   }
# 
# }



rm(list = ls())
source("ExpressionPlotsFunctions.R")
#TCGA Expression plots
#Datasets and directories and variables
datasets = c("ESCA","HNSC","LUSC","BLCA","LIHC","STAD","LGG","COAD","PAAD","READ","SKCM","LUAD")
datasets = c("BRCA")
input.dir = "~/Documents/PhD/GenderAnalysis/TCGA/Analysis/TCGAExpressionExplorerOutput/"
output.dir = "~/Documents/PhD/GenderAnalysis/TCGA/Analysis/TCGAExpressionExplorerOutput/"
comp.gene.muts = c("TP53")
gene.muts.order = list(list('Normal','Wt','Mt'))
#gene.muts.order = list(list('Wt','Mt'))

h.gene.mut.order = hash(keys=comp.gene.muts, values = gene.muts.order)
comp.other.vars = c("SAMPLE_TYPE")
#comp.other.vars = NULL
other.vars.order = list(list("Normal",'Tumor'))
#other.vars.order = list(list('Tumor'))

h.other.var.order = hash(keys=comp.other.vars, values = other.vars.order)
comp.clinical.vars = c('GENDER')
signature.file =  "~/Documents/PhD/GenderAnalysis/TP53_interactions/TP53.list.3.txt"
####Read Annotation file####
gene.annot = read.csv("~/Documents/PhD/GenderAnalysis/TCGA/Analysis/annotation.all.genes.csv",as.is = T)
gene.annot = gene.annot[,-1]
gene.annot = gene.annot[gene.annot$Chrm!='',]
gene.annot =gene.annot[order(gene.annot$Chrm,gene.annot$Entrez.Gene.ID),]
#########################

###Mutations####
print('Reading Mutations')
mutations = read.csv("~/Documents/PhD/GenderAnalysis/TCGA/Analysis/full.reduced.all.raw.TCGA.curated.mutations.csv", as.is = T)
mutations$PATIENT_ID = gsub("TCGA\\-([[:alnum:]]{2})\\-([[:alnum:]]{4}).*","\\2", mutations$SAMPLE_ID)
###############

###Clinical Data####
print('Reading Clinical')
demo = read.csv("~/Documents/PhD/GenderAnalysis/TCGA/Analysis/all.TCGA.curated.clinical.csv", as.is = T)
demo$PATIENT_ID = toupper(demo$PATIENT_ID)

###############
###Read Signatures####
sig.lists = read.delim(signature.file, sep = " ")
#####################
dir.create(paste(output.dir,sep = ""),showWarnings = F)
###Big Matrix########
all.norm.counts = NULL

###Total Results directory####
dir.create(paste(output.dir,'AllDataSets/',sep = ""),showWarnings = F)

###Start Main Loop###
print('Main Loop')
for (i in datasets)
{
  print(i)
  ###Create directories####
  dir.create(paste(output.dir,i,sep = ""),showWarnings = F)
  dir.create(paste(output.dir,i,"/ExpressionPlots",sep = ""),showWarnings = F)
  #########################
  
  ###Read Norm counts#####
  norm.counts = read.csv(paste(input.dir,i,'/Normalisation/Norm.Log.Counts.csv',sep = ""), as.is = T)
  rownames(norm.counts) = norm.counts[,1]
  norm.counts = norm.counts[,-1]
  #########################
  
  ####Prepare annotation#####
  annot = as.data.frame(cbind(SAMPLE_EXP=colnames(norm.counts)))
  annot$PATIENT_ID = gsub("TCGA\\.([[:alnum:]]{2})\\.([[:alnum:]]{4}).*","\\2", annot[,1])
  annot$SAMPLE_ID = gsub("TCGA\\.([[:alnum:]]{2})\\.([[:alnum:]]{4})\\.([[:alnum:]]{3}).*",
                         "TCGA\\-\\1\\-\\2\\-\\3", annot[,1])
  
  #Add Clinical Variables
  for (j in comp.clinical.vars){
    if (j %in% colnames(demo)){
      cols = colnames(annot)
      annot <- cbind(annot,demo[match(annot$PATIENT_ID, demo$PATIENT_ID),j]) 
      colnames(annot) = c(cols,j)
    }
  } 
  
  #Add Sample_Type#
  annot$SAMPLE_CODE = gsub("TCGA\\.([[:alnum:]]{2})\\.([[:alnum:]]{4}\\.([[:alnum:]]{2})).*","\\3", annot[,1])
  annot$SAMPLE_TYPE <- ifelse(as.integer(annot$SAMPLE_CODE)>=20, 'Control',
                              ifelse(as.integer(annot$SAMPLE_CODE)>=10, 'Normal', 'Tumor'))
  
  #Add Mutants
  filter.out = c("Silent",'Intron','IGR',"In_Frame_Ins" ,"In_Frame_Del", "lincRNA" )
  mutations %>% 
    filter(CANCER_TYPE == i,!(VARIANT_CLASSIFICATION%in%filter.out)) -> muts.dataset
  for (g in comp.gene.muts){
    g2 = paste(g,'_STATUS',sep = '')
    muts.gene = muts.dataset[which(muts.dataset$HUGO_SYMBOL==g),]
    annot[,g2] = NA
    annot[which(annot$SAMPLE_TYPE == 'Normal'), g2] = 'Normal'
    annot[which(annot$PATIENT_ID %in% muts.gene$PATIENT_ID&annot$SAMPLE_TYPE=='Tumor'),g2] = 'Mt'
    annot[which(annot$SAMPLE_TYPE=='Tumor'&is.na(annot[,g2])),g2] = 'Wt'
  }
  ################################
  
  ###Add Signatures and Genes#####
  for (s in rownames(sig.lists)){
    #Add signature
    genes = unlist(strsplit(as.character(sig.lists[s,'GENES']),',')) 
    e = gene.annot[which(gene.annot$Approved.Symbol%in%genes),'Entrez.Gene.ID']
    e.data = t(norm.counts[which(rownames(norm.counts)%in%e),])
    #a=scale(data.matrix(rowSums(e.data)))
    a=data.matrix(rowSums(e.data))
    #med.col = median(a)
    #a = scale((a)-(med.col))
    annot[,as.character(sig.lists[s,'SIGNATURE'])] = a[match(annot$SAMPLE_EXP,rownames(a)),]
    #Add individual genes
    a=data.matrix(e.data)
    colnames(a) = gene.annot[match(colnames(a),gene.annot$Entrez.Gene.ID),"Approved.Symbol"]
    for (gs in genes){
      if (gs%in%colnames(a)&&!gs%in%annot){
        med.gs = median(a[,gs])
        #annot[,as.character(gs)] = as.vector(scale(a[,gs]-(med.gs)))  
        #annot[,as.character(gs)] = as.vector(scale(a[,gs]))  
        annot[,as.character(gs)] = as.vector(a[,gs])
      }
      if(!(gs%in%colnames(a))){
        annot[,as.character(gs)] = 0
      }
      
    }
    rm(genes,e,e.data,a,med.col,med.gs)
  }
  
  ##Filter Samples
  #annot %>% filter(SAMPLE_TYPE=="Tumor") -> annot
  
  
  ################################
  #HACK FOR ONLY ONE CLINICAL VALUE (ex ONLY FEMALES)
  #annot = annot[which(annot$GENDER == 'female'),]
  ##########PLOT##################
  DoExpressionPlots(annot = annot,i = i,comp.clinical.vars = comp.clinical.vars,
                      comp.other.vars = comp.other.vars,comp.gene.muts = comp.gene.muts,
                    h.gene.mut.order = h.gene.mut.order,sig.lists = sig.lists, 
                    h.other.var.order = h.other.var.order, post.fix = "_NvWtvMt")
  ################################
  
  ###Add to Total Matrix#####
  
}


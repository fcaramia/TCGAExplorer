rm(list = ls())
####Curate Uncompressed expression file from FireHose########
source("ExpressionPlotsFunctions.R")
#TCGA Expression
#Read Files
datasets = c("ESCA","HNSC","LUSC","BLCA","LIHC","STAD","LGG","COAD","PAAD","READ","SKCM","LUAD")
datasets = c("LUAD")
input.dir = "~/Documents/PhD/Data/TCGA_2016_01_28_BROAD/Expression.Data/"
output.dir = "~/Documents/Work/Data/"
mds.contrast = c('GENDER','SAMPLE_TYPE',"TP53_STATUS")
#################
###Clinical Data####
demo = read.csv("~/Documents/PhD/GenderAnalysis/TCGA/Analysis/all.TCGA.curated.clinical.csv", as.is = T)
demo$PATIENT_ID = toupper(demo$PATIENT_ID)
DNA.MUTS = read.csv("~/Documents/PhD/GenderAnalysis/TCGA/Analysis/full.reduced.all.raw.TCGA.curated.mutations.csv", as.is=T)
DNA.MUTS %>% filter(HUGO_SYMBOL == 'TP53') -> tp53.mutants
tp53.mutants = unique(tp53.mutants$PATIENT_ID)
demo$TP53_STATUS = ifelse(demo$PATIENT_ID%in%tp53.mutants,'Mt','Wt')
###############
###Build big matrix with all samples###
do.all = F
###Big Matrix########
all.norm.counts = NULL
all.annot = NULL
all.raw.counts = NULL
###Total Results directory####
dir.create(paste(output.dir,'AllDataSets/',sep = ""),showWarnings = F)
###Housekeeping Genes#####
hk.genes = c("C1orf43","CHMP2A","EMC7","GPI","PSMB2","PSMB4","RAB7A","REEP5","SNRPD3","VCP","VPS29")
p53.reg = c('DKC1','FOXP3','HDAC8','LAS1L','TAF1','YY2','OGT','DDX3X','AMER1','BMX','CUL4B','DDX53','HUWE1','MAGEA2','MCTS1','NOX1','PSMD10','UBE2A','UTP14A','UXT','XIAP','AIFM1','G6PD','CD40LG','IL2RG','SH2D1A','TLR8','AR')
####Read Annotation file####
gene.annot = read.csv("~/Documents/PhD/GenderAnalysis/TCGA/Analysis/annotation.all.genes.csv",as.is = T)
gene.annot = gene.annot[,-1]
gene.annot = gene.annot[gene.annot$Chrm!='',]
gene.annot =gene.annot[order(gene.annot$Chrm,gene.annot$Entrez.Gene.ID),]
#########################
for (i in datasets)
{
  print(i)
  ###Create directories####
  dir.create(paste(output.dir,i,sep = ""),showWarnings = F)
  dir.create(paste(output.dir,i,"/Normalisation",sep = ""),showWarnings = F)
  #########################
  
  ###Read Raw counts#####
  raw.counts = read.csv(paste(input.dir,i,'/RawExpression.csv',sep = ""))
  rownames(raw.counts) = raw.counts[,1]
  raw.counts = raw.counts[,-1]
  #########################
  #Keep ENTREZ.ID
  rownames(raw.counts) = gsub(".*\\|(.*)","\\1", rownames(raw.counts))
  #########################
  raw.counts = data.matrix(raw.counts)
  ###Normalise#############
  logCPM = DoCPMNorm(raw.counts)
  #########################
  ##Save Norm counts####
  write.csv(logCPM,paste(sep="",output.dir,i,"/Normalisation","/Norm.Log.Counts.csv"))
  ##############
  ###Compute Annotations####
  annot = as.data.frame(cbind(SAMPLE_ID=colnames(raw.counts)))
  annot$PATIENT_ID = gsub("TCGA\\.([[:alnum:]]{2})\\.([[:alnum:]]{4}).*","\\2", annot[,1])
  for (j in mds.contrast){
    if (j %in% colnames(demo)){
      cols = colnames(annot)
      annot <- cbind(annot,demo[match(annot$PATIENT_ID, demo$PATIENT_ID),j]) 
      colnames(annot) = c(cols,j)
    }
  } 
  annot$DATA_SET = i
  annot$SAMPLE_CODE = gsub("TCGA\\.([[:alnum:]]{2})\\.([[:alnum:]]{4})\\.([[:alnum:]]{2}).*","\\3", annot[,1])
  annot$SAMPLE_TYPE <- ifelse(as.integer(annot$SAMPLE_CODE)>=20, 'Control',
                ifelse(as.integer(annot$SAMPLE_CODE)>=10, 'Normal', 'Tumor')
  )
  ##########################
  ###Do Plots####
  # DoNormPlots(raw.data = raw.counts, norm.data = logCPM,
  #             dir = paste(output.dir,i,"/Normalisation", sep ="" ), tittle = i,
  #             annot = annot, color.cols = mds.contrast)
  # DoMultMDSPlots(norm.data = logCPM, dir = paste(output.dir,i,"/Normalisation", sep ="" ), tittle = i,
  #                annot = annot, color.cols = mds.contrast )
  # ###############

  if (do.all == T){
    if(is.null(all.raw.counts)){
    #  all.norm.counts=logCPM
      all.annot = annot
      all.raw.counts = raw.counts
    }else{
     # genes.in.common = intersect(rownames(all.norm.counts),rownames(logCPM))
      #all.norm.counts <- cbind(all.norm.counts[genes.in.common , , drop=FALSE], logCPM[genes.in.common , , drop=FALSE])
      all.annot <- rbind(all.annot, annot)
      genes.in.common = intersect(rownames(all.raw.counts),rownames(raw.counts))
      all.raw.counts <- cbind(all.raw.counts[genes.in.common , , drop=FALSE], raw.counts[genes.in.common , , drop=FALSE])
    }
  }
}
if (do.all == T){
  #ind.samples = c(sample(1:n,500,replace = F))
  #ind.genes = as.character(gene.annot[which(gene.annot$Approved.Symbol%in%hk.genes),'Entrez.Gene.ID'])
  #ind.genes = c(sample(1:15000,50,replace = F))
  #samples = all.annot[ind.samples,'SAMPLE_ID']
  #samples = as.vector(all.annot[which(all.annot$SAMPLE_TYPE=='Normal'),'SAMPLE_ID'])
  
  
  
  
  aux <- factor(paste(all.annot$DATA_SET,all.annot$SAMPLE_TYPE,sep="."))
  design <- model.matrix(~0+aux)
  colnames(design) <- levels(aux)
  
  source("ExpressionPlotsFunctions.R")
  #all.norm.counts <- NormaliseVoom(all.raw.counts,design = design)
  all.norm.counts = DoCPMNorm(all.raw.counts)
  n = ncol(all.norm.counts)
  
  #sub.mat = all.norm.counts[ind.genes,samples]
  #DoMultMDSPlots(norm.data = sub.mat, dir = paste(output.dir,"AllDataSets", sep ="" ), tittle = 'All_DATA_SETS',
  #               annot = all.annot[which(all.annot$SAMPLE_ID%in%samples),], color.cols = c("DATA_SET") )
  
  write.csv(all.norm.counts,paste(sep="",output.dir,"AllDataSets","/Norm.Log.Counts.csv"))
}

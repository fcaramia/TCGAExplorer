rm(list = ls())
library(data.table)
library(openxlsx)
library(dplyr)
library(genefu)
source("~/Documents/Work/Projects/Sherene/Combined_Plotting/survival_functions.R")
source('PropensityScoresFunctions.R')
source("ExpressionPlotsFunctions.R")

input.dir = "~/Documents/PhD/GenderAnalysis/TCGA/Analysis/TCGAExpressionExplorerOutput/"
output_dir="~/Documents/PhD/GenderAnalysis/TCGA/Analysis/TCGAExpressionExplorerOutput/SurvivalAnalysis2/"
if (!file.exists(output_dir)){
  dir.create(file.path(output_dir))
}


adjust.hr = c('age','race','PATHOLOGIC_STAGE','SMOKING_STATUS')
adjust.hr = NULL

expression = fread("~/Documents/PhD/Data/TCGA_Xena/Expression/EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena")
datasets = c("ESCA","HNSC","LUSC","BLCA","LIHC","STAD","LGG","COAD","PAAD","READ","SKCM",'LUAD')
clinical = read.xlsx("~/Documents/PhD/Data/TCGA_CLINICAL/mmc1.xlsx", sheet = 1)
clinical = filter(clinical, type%in%datasets)
clinical$PATIENT_ID = gsub("TCGA\\-([[:alnum:]]{2})\\-([[:alnum:]]{4}).*","\\2",clinical$bcr_patient_barcode)
clinical <- clinical %>% mutate( age = ifelse(age_at_initial_pathologic_diagnosis >=45, "old", "young"))
demo = read.csv("~/Documents/PhD/GenderAnalysis/TCGA/Analysis/all.TCGA.curated.clinical.csv", as.is = T)
demo$PATIENT_ID = toupper(demo$PATIENT_ID)
demo= filter(demo, CANCER_TYPE%in%datasets)
clinical = full_join(clinical, demo)


####Read Annotation file####
gene.annot = read.csv("~/Documents/PhD/GenderAnalysis/TCGA/Analysis/annotation.all.genes.csv",as.is = T)
gene.annot = gene.annot[,-1]
gene.annot = gene.annot[gene.annot$Chrm!='',]
gene.annot = gene.annot[order(gene.annot$Chrm,gene.annot$Entrez.Gene.ID),]
#########################



##P53 Interactors
p53.int = read.csv("~/Documents/PhD/GenderAnalysis/TP53_interactions/Candidates.2019.02.01.csv")
p53.low = as.vector(p53.int[which(p53.int$experimentally_determined_interaction>=0.3),'GeneSymbol'])
p53.low = as.vector(p53.int[,'GeneSymbol'])

p53.int$EntrezID = gene.annot[match(p53.int$GeneSymbol,gene.annot$Approved.Symbol),'Entrez.Gene.ID']



#####Filter p53 genes
expression <- subset(expression, select=which(!duplicated(names(expression)))) 
expression = filter(expression,sample%in%p53.int$GeneSymbol)
rownames(expression) = expression$sample
expression$sample = NULL
expression = as.data.frame(t(expression))
expression$Sample = rownames(expression)
expression$SampleType = gsub("TCGA\\-([[:alnum:]]{2})\\-([[:alnum:]]{4})\\-([[:alnum:]]{2})","\\3",rownames(expression))
expression$PATIENT_ID = gsub("TCGA\\-([[:alnum:]]{2})\\-([[:alnum:]]{4}).*","\\2",rownames(expression))
expression = filter(expression,PATIENT_ID%in%clinical$PATIENT_ID,SampleType!='11')
mat = full_join(expression,clinical)

###Compute signature 

sig = p53.int[which(p53.int$GeneSymbol%in%colnames(mat)&p53.int$GeneSymbol%in%p53.low),c("GeneSymbol",'EntrezID')]
colnames(sig) = c('probe','EntrezGene.ID')
sig$coefficient = 1

score = sig.score(x = sig,data = mat,annot = sig)
mat$p53.string.sig = score$score

DNA.muts = fread("~/Documents/PhD/GenderAnalysis/TCGA/Analysis/full.reduced.all.raw.TCGA.curated.mutations.csv")
DNA.muts %>% filter(CANCER_TYPE%in%datasets,HUGO_SYMBOL=='TP53') -> DNA.muts

filter.out = c("Silent",'Intron','IGR',"In_Frame_Ins" ,"In_Frame_Del", "lincRNA" )

filter(DNA.muts,!(VARIANT_CLASSIFICATION%in%filter.out)) -> DNA.muts
p53.mutants = unique(DNA.muts$PATIENT_ID)

mat$TP53.STATUS = ifelse(mat$PATIENT_ID%in%p53.mutants,'Mt-p53','Wt-p53')


####Create MAtrix for individual cancer.

expression.list = list()
for (cancer in datasets){
  print(cancer)
  norm.counts = data.matrix(fread(paste(input.dir,cancer,'/Normalisation/Norm.Log.Counts.csv',sep = "")))
  rownames(norm.counts) = norm.counts[,1]
  norm.counts = norm.counts[,-1]
  #Use counts with clinical data
  annot = as.data.frame(cbind(SAMPLE_EXP=colnames(norm.counts)))
  annot$PATIENT_ID = gsub("TCGA\\.([[:alnum:]]{2})\\.([[:alnum:]]{4}).*","\\2", annot[,1])
  annot$SAMPLE_ID = gsub("TCGA\\.([[:alnum:]]{2})\\.([[:alnum:]]{4})\\.([[:alnum:]]{3}).*",
                         "TCGA\\-\\1\\-\\2\\-\\3", annot[,1])
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
  
  
  annot$TP53.STATUS = NA
  annot[which(annot$SAMPLE_TYPE == 'Normal'), 'TP53.STATUS'] = 'Normal'
  annot[which(annot$PATIENT_ID %in% p53.mutants&annot$SAMPLE_TYPE=='Tumor'),'TP53.STATUS'] = 'Mt'
  annot[which(annot$SAMPLE_TYPE=='Tumor'&is.na(annot[,'TP53.STATUS'])),'TP53.STATUS'] = 'Wt'

  ####ENTREZ.ID to SYMBOL
  norm.counts = norm.counts[which(rownames(norm.counts)%in%p53.int$EntrezID),]
  rownames(norm.counts) = gene.annot[match(rownames(norm.counts),gene.annot[,'Entrez.Gene.ID']),'Approved.Symbol']
  sig.exp = data.frame(t(norm.counts))
  
  sig = p53.int[which(p53.int$GeneSymbol%in%colnames(sig.exp)&p53.int$GeneSymbol%in%p53.low),c("GeneSymbol",'EntrezID')]
  colnames(sig) = c('probe','EntrezGene.ID')
  sig$coefficient = 1
  
  score = sig.score(x = sig,data = sig.exp,annot = sig)
  sig.exp$p53.string.sig = score$score
  
  sig.exp$PATIENT_ID = gsub("TCGA\\.([[:alnum:]]{2})\\.([[:alnum:]]{4}).*","\\2", rownames(sig.exp))
  
  sig.exp = full_join(sig.exp,annot)
  ###Only tumors
  sig.exp = filter(sig.exp,SAMPLE_TYPE=='Tumor')
  cancer.mat = left_join(sig.exp,clinical)

  expression.list[[cancer]] = cancer.mat
}




####Plot densities#####
#DoDensityPlotGenderRMAF(dat = mat, dat.col = 'XIST',fill.col = 'gender',title.txt = 'XIST Distribution')


########################


# ####Discretize p53 string mutations by patient
# DNA.muts = fread("~/Documents/PhD/GenderAnalysis/TCGA/Analysis/full.reduced.all.raw.TCGA.curated.mutations.csv")
# DNA.muts.p53.string = filter(DNA.muts,HUGO_SYMBOL%in%as.vector(sig$probe ))
# as.data.frame(DNA.muts.p53.string %>% group_by(PATIENT_ID) %>% dplyr::summarise(N.P = n())) -> desc.p53.string.muts.numbers
# mat$p53.string.muts.desc = desc.p53.string.muts.numbers[match(mat$PatientID,desc.p53.string.muts.numbers$PATIENT_ID),'N.P']
# mat[which(is.na(mat$p53.string.muts.desc)),'p53.string.muts.desc'] = 0
# 
# mat$p53.string.muts.desc.d = ifelse(mat$p53.string.muts.desc==0,'Wt',ifelse(mat$p53.string.muts.desc>6,NA,ifelse(mat$p53.string.muts.desc>2,'3 to 6','1 to 2')))
# 
# pdf(paste(output_dir,"/TCGA_p53_set_mutations_survival_os.pdf",sep=""))
# 
# mat$comb = paste(mat$p53.string.muts.desc.d,mat$p53.status,sep = '.')
# mat2 = mat[!is.na(mat$p53.string.muts.desc.d)&mat$p53.string.muts.desc.d!='Wt',]
# 
# survplot(mat2,quantcol="p53.string.muts.desc.d",survcol="OS",timecol="OS.time",subset.mat = NULL, show.conf.int = F,
#          contcol = 'p53.string.muts.desc.d',title= paste("TCGA p53 string muts"),percentile.colors=c("blue","grey","red",'black'),
#          time.mult=365,use.cont = F,ylab="Overall Survival",max.x = 5)
# 
# 
# survplot(mat2,quantcol="comb",survcol="OS",timecol="OS.time",subset.mat = NULL, show.conf.int = F,
#          contcol = 'comb',title= paste("TCGA Gender"),percentile.colors=c("blue","grey","red",'black','yellow','green','orange','purple'),
#          time.mult=365,use.cont = F,ylab="Overall Survival",max.x = 5)
# 
# 
# dev.off()
# 
# 
# 
# ####Do plots
# 
# 
# genes = c(as.vector(sig$probe),'p53.string.sig')
# #genes = c('p53.string.sig')
# 
# ###Male and female 
# 
# 
# tiff(paste(output_dir,"/TCGA_42_mts_p53_survival_dss.tiff",sep=""),width = 2500,height = 2800,res = 500)
# survplot(mat,quantcol="p53.status",survcol="DSS",timecol="DSS.time",subset.mat = NULL, show.conf.int = F,surv.line = NULL,
#          contcol = 'p53.status',title= paste("TCGA Gender"),percentile.colors=c("turquoise1",'darkorchid1'), mark.size = 0,
#          time.mult=33.3,use.cont = F,ylab="Disease Specific Survival",max.x = 42, line.type = c(1,1))
# dev.off()
# mat$comb = paste(mat$p53.status,mat$gender,sep = '.')
# 
# tiff(paste(output_dir,"/TCGA_5_yrs_p53_gender_survival_dss.tiff",sep=""),width = 2500,height = 2800,res = 500)
# survplot(mat,quantcol="comb",survcol="DSS",timecol="DSS.time",subset.mat = NULL, show.conf.int = F,
#          contcol = 'comb',title= paste("TCGA"),percentile.colors=c('red','black','lightcoral','mediumspringgreen'),
#          time.mult=365,use.cont = F,ylab="Disease Specific Survival",max.x = 5, line.type = c(2,2,3,1))
# 
# dev.off()
# 
# #####Only Females and males
# mat2 = mat[mat$gender=='FEMALE',]
# mat2$comb = paste(mat2$p53.status,mat2$gender,sep = '.')
# 
# tiff(paste(output_dir,"/TCGA_14_yrs_p53_female_survival_dss.tiff",sep=""),width = 2000,height = 2000,res = 270)
# survplot(mat2,quantcol="comb",survcol="DSS",timecol="DSS.time",subset.mat = NULL, show.conf.int = F,
#          contcol = 'p53.status',title= paste("TCGA"),percentile.colors=c('red','black'),
#          time.mult=365,use.cont = F,ylab="Disease Specific Survival",max.x = 14, line.type = c(3,3))
# 
# dev.off()
# 
# 
# mat2 = mat[mat$gender=='MALE',]
# mat2$comb = paste(mat2$p53.status,mat2$gender,sep = '.')
# 
# tiff(paste(output_dir,"/TCGA_5_yrs_gender_p53_survival_os.tiff",sep=""),width = 2500,height = 2800,res = 400)
# survplot(mat2,quantcol="comb",survcol="OS",timecol="OS.time",subset.mat = NULL, show.conf.int = F,
#          contcol = 'p53.status',title= paste("TCGA"),percentile.colors=c('red','black'),
#          time.mult=365,use.cont = F,ylab="Overall Survival",max.x = 5, line.type = c(1,1))
# 
# dev.off()
# 
# 
# 
# # for (g in genes){
# # 
# #   gene.d = paste(g,'d',sep = '.')
# #   mat[,gene.d] = ifelse(mat[,g]>=median(mat[,g],na.rm = T),paste('hi',g,sep = '-'),paste('lo',g,sep = '-'))
# #   mat$comb = paste(mat$gender,mat[,gene.d],sep = '.')
# # 
# #   pdf(paste(output_dir,"/TCGA_",g,"_gender_survival_os.pdf",sep=""))
# #   mat2 = mat[which(!is.na(mat[,g])),]
# #   survplot(mat2,quantcol="comb",survcol="OS",timecol="OS.time",subset.mat = NULL, show.conf.int = F,
# #            contcol = g,title= paste("TCGA ",g),percentile.colors=c("blue","grey","red",'black'),
# #            time.mult=365,use.cont = F,ylab="Overall Survival",max.x = 15)
# # 
# #   dev.off()
# # 
# # 
# # }
# 
# #####p53 Wt vs p53 Mt
# ###For Sue
# ####p53 signature only
# bckup = mat
# g= "XIST"
# demo = 'DSS'
# years = 14
# y.lab = ifelse(demo=='DSS','Disease Specific Survival',ifelse(demo=="OS","Overall Survival","Disease Free Survival"))
# evantcolname = ifelse(demo=='DSS','DSS',ifelse(demo=="OS","OS","e.dfs"))
# timecolname = ifelse(demo=='DSS','DSS.time',ifelse(demo=="OS","OS.time","t.dfs"))
# 
# 
# gene.d = paste(g,'d',sep = '.')
# mat = mat[which(!is.na(mat[,g])&mat$gender=='FEMALE'),]
# 
# mat[,gene.d] = ifelse(mat[,g]>=median(mat[,g],na.rm = T),paste('hi',g,sep = '-'),paste('lo',g,sep = '-'))
# mat$comb = paste(mat$p53.status,mat[,gene.d],mat$gender,sep = '.')
# 
# tiff(paste(output_dir,"/TCGA_",years,"_yrs__p53_status_XIST_female_survival_",demo,".tiff",sep=""),width = 2500,height = 2800,res = 500)
# survplot(mat,quantcol="comb",survcol=evantcolname,timecol=timecolname,subset.mat = NULL, show.conf.int = F, line.type = c(2,1,2,1),
#          contcol = g,title= paste(g),percentile.colors=c('blue','black','magenta','red'), mark.size = 0.5, 
#          line.width = 2,
#          time.mult=365,use.cont = F,ylab=y.lab,max.x = years,print.significant.only = F)
# dev.off()
# mat = bckup



#Do by cancer type
mat2 = mat
for (c in c(datasets,'All')){
  print(c)
  adjust.tmp = NULL
  if (c == 'All'){
    mat = mat2
    #adjust.tmp = c('CANCER_TYPE')
  }else{
    mat = expression.list[[c]]  
  }

  for (cf in adjust.hr){
    na.cont1 = length(mat[which(mat[,cf] == "[Not Available]"),cf])
    na.cont1 = length(mat[which(is.na(mat[,cf])),cf]) + na.cont1
    na.cont1 = length(mat[which(mat[,cf] == "[Not Evaluated]"),cf])  + na.cont1
    na.cont1 = length(mat[which(mat[,cf] == "[Unknown]"),cf]) + na.cont1
    half = dim(mat)[1]/4
    
    
    if(na.cont1 < half){
      adjust.tmp = c(adjust.tmp,cf)
    }
    
  }
  
  c.dir = paste(output_dir,c,sep = '')
  if (!file.exists(c.dir)){
    dir.create(file.path(c.dir))
  }
  pdf(paste(c.dir,"/TCGA_",c,"_p53_sig_survival.pdf",sep=""))
  
  for( d in c('OS','DSS') ){
  
    for (t in c(5,14)){
      bckup = mat
      g= "p53.string.sig"
      demo = d
      years = t
      y.lab = ifelse(demo=='DSS','Disease Specific Survival',ifelse(demo=="OS","Overall Survival","Disease Free Survival"))
      evantcolname = ifelse(demo=='DSS','DSS',ifelse(demo=="OS","OS","e.dfs"))
      timecolname = ifelse(demo=='DSS','DSS.time',ifelse(demo=="OS","OS.time","t.dfs"))
      gene.d = paste(g,'d',sep = '.')
      
      # chisq.test( table( mat$gender, mat$TP53.STATUS ) )
      # chisq.test( table( mat$gender, mat$PATHOLOGIC_STAGE ) )
      # 
      # chisq.test( table( mat$TP53.STATUS, mat$PATHOLOGIC_STAGE ) )
      
      
      # mat = mat[which(!is.na(mat[,g])),]
      # mat[,gene.d] = ifelse(mat[,g]>=median(mat[,g],na.rm = T),paste('hi',g,sep = '-'),paste('lo',g,sep = '-'))
      # mat$Var = as.factor(paste(mat$gender,sep = '.'))
      # 
      # 
      # survplot(mat,quantcol='gender',survcol=evantcolname,timecol=timecolname,subset.mat = NULL, show.conf.int = F, line.type = c(2,1,2,1),
      #          contcol = 'gender',title= paste(g),percentile.colors=c('blue','black','magenta','red'), mark.size = 0.5, 
      #          line.width = 2,
      #          time.mult=365,use.cont = F,ylab=y.lab,max.x = years,print.significant.only = F , adjust.hr.by = adjust.tmp)
      # 
      # mat = bckup
      # 
      # mat = mat[which(!is.na(mat[,g])),]
      # mat[,gene.d] = ifelse(mat[,g]>=median(mat[,g],na.rm = T),paste('hi',g,sep = '-'),paste('lo',g,sep = '-'))
      # mat$Var = as.factor(paste(mat$TP53.STATUS,sep = '.'))
      # 
      # 
      # survplot(mat,quantcol='TP53.STATUS',survcol=evantcolname,timecol=timecolname,subset.mat = NULL, show.conf.int = F, line.type = c(2,1,2,1),
      #          contcol = 'TP53.STATUS',title= paste(g),percentile.colors=c('blue','black','magenta','red'), mark.size = 0.5, 
      #          line.width = 2,
      #          time.mult=365,use.cont = F,ylab=y.lab,max.x = years,print.significant.only = F , adjust.hr.by = adjust.tmp)
      # 
      # mat = bckup
      # 
      # mat = mat[which(!is.na(mat[,g])),]
      # mat[,gene.d] = ifelse(mat[,g]>=median(mat[,g],na.rm = T),paste('hi',g,sep = '-'),paste('lo',g,sep = '-'))
      # mat$Var = as.factor(paste(mat$TP53.STATUS,mat$gender,sep = '.'))
      # 
      # 
      # survplot(mat,quantcol='Var',survcol=evantcolname,timecol=timecolname,subset.mat = NULL, show.conf.int = F, line.type = c(2,1,2,1),
      #          contcol = 'Var',title= paste(g),percentile.colors=c('blue','black','magenta','red'), mark.size = 0.5, 
      #          line.width = 2,
      #          time.mult=365,use.cont = F,ylab=y.lab,max.x = years,print.significant.only = F , adjust.hr.by = adjust.tmp)
      # 
      # mat = bckup
      # 
      # 
      # 
      # mat = mat[which(!is.na(mat[,g])),]
      # mat[,gene.d] = ifelse(mat[,g]>=median(mat[,g],na.rm = T),paste('hi',g,sep = '-'),paste('lo',g,sep = '-'))
      # mat$Var = as.factor(paste(mat[,gene.d],sep = '.'))
      # 
      # 
      # survplot(mat,quantcol="Var",survcol=evantcolname,timecol=timecolname,subset.mat = NULL, show.conf.int = F, line.type = c(2,1,2,1),
      #          contcol = g,title= paste(g),percentile.colors=c('blue','black','magenta','red'), mark.size = 0.5, 
      #          line.width = 2,
      #          time.mult=365,use.cont = F,ylab=y.lab,max.x = years,print.significant.only = F , adjust.hr.by = adjust.tmp)
      # 
      # mat = bckup
      # 
      # 
      # mat = mat[which(!is.na(mat[,g])),]
      # mat[,gene.d] = ifelse(mat[,g]>=median(mat[,g],na.rm = T),paste('hi',g,sep = '-'),paste('lo',g,sep = '-'))
      # mat$Var = as.factor(paste(mat$TP53.STATUS,mat[,gene.d],sep = '.'))
      # 
      # 
      # survplot(mat,quantcol="Var",survcol=evantcolname,timecol=timecolname,subset.mat = NULL, show.conf.int = F, line.type = c(2,1,2,1),
      #          contcol = "Var",title= paste(g),percentile.colors=c('blue','black','magenta','red'), mark.size = 0.5, 
      #          line.width = 2,
      #          time.mult=365,use.cont = F,ylab=y.lab,max.x = years,print.significant.only = F , adjust.hr.by = adjust.tmp)
      # 
      # mat = bckup
      # 
      # 
      # mat = mat[which(!is.na(mat[,g])),]
      # mat[,gene.d] = ifelse(mat[,g]>=median(mat[,g],na.rm = T),paste('hi',g,sep = '-'),paste('lo',g,sep = '-'))
      # mat$Var = as.factor(paste(mat[,gene.d],mat$gender,sep = '.'))
      # 
      # 
      # survplot(mat,quantcol="Var",survcol=evantcolname,timecol=timecolname,subset.mat = NULL, show.conf.int = F, line.type = c(2,1,2,1),
      #          contcol = "Var",title= paste(g),percentile.colors=c('blue','black','magenta','red'), mark.size = 0.5, 
      #          line.width = 2,
      #          time.mult=365,use.cont = F,ylab=y.lab,max.x = years,print.significant.only = F, adjust.hr.by = adjust.tmp)
      # 
      # mat = bckup
      # 
      # mat = mat[which(!is.na(mat[,g])&mat$gender=='FEMALE'),]
      # mat[,gene.d] = ifelse(mat[,g]>=median(mat[,g],na.rm = T),paste('hi',g,sep = '-'),paste('lo',g,sep = '-'))
      # mat$Var = as.factor(paste(mat$TP53.STATUS,mat[,gene.d],mat$gender,sep = '.'))
      # 
      # 
      # survplot(mat,quantcol="Var",survcol=evantcolname,timecol=timecolname,subset.mat = NULL, show.conf.int = F, line.type = c(2,1,2,1),
      #          contcol = "Var",title= paste(g),percentile.colors=c('blue','black','magenta','red'), mark.size = 0.5, 
      #          line.width = 2,
      #          time.mult=365,use.cont = F,ylab=y.lab,max.x = years,print.significant.only = F, adjust.hr.by = adjust.tmp)
      # 
      # mat = bckup
      # 
      # mat = mat[which(!is.na(mat[,g])&mat$gender=='MALE'),]
      # mat[,gene.d] = ifelse(mat[,g]>=median(mat[,g],na.rm = T),paste('hi',g,sep = '-'),paste('lo',g,sep = '-'))
      # mat$Var = as.factor(paste(mat$TP53.STATUS,mat[,gene.d],mat$gender,sep = '.'))
      # 
      # 
      # survplot(mat,quantcol="Var",survcol=evantcolname,timecol=timecolname,subset.mat = NULL, show.conf.int = F, line.type = c(2,1,2,1),
      #          contcol = "Var",title= paste(g),percentile.colors=c('blue','black','magenta','red'), mark.size = 0.5, 
      #          line.width = 2,
      #          time.mult=365,use.cont = F,ylab=y.lab,max.x = years,print.significant.only = F, adjust.hr.by = adjust.tmp)
      # 
      # mat = bckup
    }  
    
  }
  dev.off()
}




  
  
####ALl genes
genes = c(as.vector(sig$probe))

pdf(paste(output_dir,"/TCGA_42_mts_females_Wt-p53_all_genes_survival_dss.pdf",sep=""))

#tiff(paste(output_dir,"/TCGA_5_yrs_Mt-p53_STATUS_FEMALE_survival_50-50prct.tiff",sep=""),width = 1200,height = 1200,res = 200)
for (g in genes){

  bckup = mat
  
  gene.d = paste(g,'d',sep = '.')
  mat[,gene.d] = ifelse(mat[,g]>=median(mat[,g],na.rm = T),paste('hi',g,sep = '-'),paste('lo',g,sep = '-'))
  mat = mat[which(!is.na(mat[,g])&mat$gender=='FEMALE'&mat$p53.status=='Wt-p53'),]
  mat$comb = paste(mat$p53.status,mat[,gene.d],mat$gender,sep = '.')
  
  
  survplot(mat,quantcol="comb",survcol="DSS",timecol="DSS.time",subset.mat = NULL, show.conf.int = F, line.type = c(3,3,3,1),
           contcol = g,title= paste(g),percentile.colors=c('blue','darkolivegreen','darkolivegreen','deepskyblue1'), mark.size = 0, 
           line.width = 2,
           time.mult=33.3,use.cont = F,ylab="Disease Specific Survival",max.x = 42,print.significant.only = F)
  
  mat = bckup
}
dev.off()
# 
# 
# demo = "OS"
# 
# y.lab = ifelse(demo=='RFS','Relapse Free Survival',ifelse(demo=="OS","Overall Survival","Disease Free Survival"))
# evantcolname = ifelse(demo=='RFS','PFI',ifelse(demo=="OS","OS","DSS"))
# timecolname = ifelse(demo=='RFS','PFI.time',ifelse(demo=="OS","OS.time","DSS.time"))
# 
# for (g in genes){
#   title = paste(g,' Survival:\nTCGA')
#   
#   
#   pdf(paste(output_dir,"/TCGA_",g,'_detailed_',demo,"_survival_percentile.pdf",sep=""),width = 5,height = 4)
#   
#   survplot_percentiles(mat,quantcol=g,survcol=evantcolname,
#                        timecol=timecolname,ylab = y.lab,use.cont = F,
#                        contcol = g,title= title,time.mult=365
#                        ,subset.mat=NULL,max.x = 15,
#                        percentiles.abs=F)
#   
#   dev.off()
#   
#   pdf(paste(output_dir,"/TCGA_",g,'_detailed_',demo,"Mt-p53_survival_percentile.pdf",sep=""),width = 5,height = 4)
#   
#   survplot_percentiles(mat[which(mat$p53.status=='Mt-p53'),],quantcol=g,survcol=evantcolname,
#                        timecol=timecolname,ylab = y.lab,use.cont = F,
#                        contcol = g,title= title,time.mult=365
#                        ,subset.mat=NULL,max.x = 15,
#                        percentiles.abs=F)
#   
#   dev.off()
#   
#   pdf(paste(output_dir,"/TCGA_",g,'_detailed_',demo,"Wt-p53_survival_percentile.pdf",sep=""),width = 5,height = 4)
#   
#   survplot_percentiles(mat[which(mat$p53.status=='Wt-p53'),],quantcol=g,survcol=evantcolname,
#                        timecol=timecolname,ylab = y.lab,use.cont = F,
#                        contcol = g,title= title,time.mult=365
#                        ,subset.mat=NULL,max.x = 15,
#                        percentiles.abs=F)
#   
#   dev.off()
#   
#   
#   
# }


#Try EM mutations carriers in multiple p53 string genes. 
RMAF = read.csv('~/Documents/PhD/GenderAnalysis/TCGA/Analysis/TCGAExpressionExplorerOutput/AllDataSets/Computed.All.RMAF.csv')
RMAF %>% filter(HUGO_SYMBOL%in%as.vector(sig$probe)) -> RMAF
#Clean low reads
RMAF$SUM = RMAF$REF_COUNT + RMAF$ALT_COUNT

RMAF %>% mutate(RMAF = ifelse(SUM>=10,ALT_COUNT/SUM,0.0)) %>% filter(!(EXPRS.GENE=='NOT.EXPRS'&SUM<10)) %>%
  filter(!(SUM<10&Z.SCORE>-2))-> RMAF

RMAF[which(RMAF$Z.SCORE<=-4&RMAF$SUM<10),'RMAF'] = 1.0
#Use mutations present in 10+ patients
RMAF %>% group_by(HUGO_SYMBOL, PATIENT_ID) %>% dplyr::summarise(s=max(RMAF)) %>% 
  dplyr::mutate(U_ID = paste(HUGO_SYMBOL, PATIENT_ID,s,sep = '.')) -> aux

RMAF$U_ID = paste(RMAF$HUGO_SYMBOL, RMAF$PATIENT_ID,RMAF$RMAF,sep = '.')
RMAF %>% filter(U_ID %in% aux$U_ID) -> RMAF

RMAF %>% group_by(HUGO_SYMBOL, PATIENT_ID) %>% dplyr::summarise(s=max(START_POSITION)) %>% 
  dplyr::mutate(U_ID = paste(HUGO_SYMBOL, PATIENT_ID,s,sep = '.')) -> aux

RMAF$U_ID = paste(RMAF$HUGO_SYMBOL, RMAF$PATIENT_ID,RMAF$START_POSITION,sep = '.')
RMAF %>% filter(U_ID %in% aux$U_ID) -> RMAF

RMAF %>% group_by(HUGO_SYMBOL) %>% dplyr::summarise(n=n()) %>% filter(n>=10) -> aux
RMAF %>% filter(HUGO_SYMBOL %in% aux$HUGO_SYMBOL) -> RMAF
rm(aux)

#RMAF$LOGIT.RMAF = logit(x = RMAF$RMAF,min = 0-0.1 , max = 1+0.1)
#RMAF[which(RMAF$Z.SCORE<=-4&RMAF$SUM<5),'RMAF'] = 1.0
#DoDensityPlotGenderRMAF(dat = RMAF,dat.col = 'RMAF',fill.col = 'GENDER',title.txt = 'RMAF')
#DoDensityPlotGenderRMAF(dat = RMAF,dat.col = 'Z.SCORE',fill.col = 'GENDER',title.txt = 'Z-Scores')
#Define Expressed Mutation RMAF
#RMAF$EM.RMAF = ifelse(RMAF$RMAF>0.3,1,0)
#RMAF$SM.RMAF = ifelse(RMAF$RMAF<=0.3,1,0)

RMAF$EM.RMAF = ifelse(RMAF$RMAF>=0.75,1,0)
RMAF$SM.RMAF = ifelse(RMAF$RMAF<=0.20,1,0)

###Filter EM muts
RMAF %>% filter(EM.RMAF==1) -> RMAF

as.data.frame(RMAF %>% group_by(PATIENT_ID) %>% dplyr::summarise(N.P = n())) -> EM.PATS
mat$p53.string.em.muts = EM.PATS[match(mat$PatientID,EM.PATS$PATIENT_ID),'N.P']
mat[which(is.na(mat$p53.string.em.muts)),'p53.string.em.muts'] = 0

mat$p53.string.EM.desc = ifelse(mat$p53.string.em.muts==0,'Wt','1+')

mat$comb = paste(mat$p53.string.EM.desc,mat$p53.status,mat$gender,sep = '.')
#mat.plot = mat[which(mat$p53.string.EM.desc!='Wt/1'),]
mat.plot = mat[mat$p53.status=='Wt-p53',]
survplot(mat.plot,quantcol="p53.string.EM.desc",survcol="DSS",timecol="DSS.time",subset.mat = NULL, show.conf.int = F,
         contcol = 'p53.string.EM.desc',title= paste("TCGA p53 string EM"),percentile.colors=c("blue","grey","red",'black'),
         time.mult=365,use.cont = F,ylab="Overall Survival",max.x = 5)

survplot(mat.plot,quantcol="comb",survcol="DSS",timecol="DSS.time",subset.mat = NULL, show.conf.int = F,
         contcol = 'comb',title= paste("TCGA Gender"),percentile.colors=c("blue","grey","red",'black','yellow','green','orange','purple'),
         time.mult=365,use.cont = F,ylab="Overall Survival",max.x = 5)




#Survival based on expressed mutations
pdf(paste(output_dir,"/TCGA_status_survival_Expressed_Muts.pdf",sep=""))
for (g in genes){
  
  gene.d = paste(g,'d',sep = '.')
  RMAF %>% filter(HUGO_SYMBOL==g) -> RMAF.g
  # mat[,gene.d] = ifelse(mat[,g]>=quantile(mat[,g],na.rm = T)[4],paste('hi',g,sep = '-'),
  #                       ifelse(mat[,g]<=quantile(mat[,g],na.rm = T)[2],paste('lo',g,sep = '-'),NA))
  #mat[,gene.d] = ifelse(mat[,g]>=median(mat[,g],na.rm = T),paste('hi',g,sep = '-'),paste('lo',g,sep = '-'))
  mat$g.mut = ifelse(mat$PatientID%in%RMAF.g$PATIENT_ID,paste('Mt',g,sep = '-'),paste("Wt",g,sep = '-'))
  mat$comb = paste(mat$p53.status,mat[,'g.mut'],sep = '.')
  
  
  #mat2 = mat[which(!is.na(mat[,g])&!is.na(mat[,gene.d])&mat[,"p53.status"]=='Wt-p53'),]
  survplot(mat,quantcol="comb",survcol="DSS",timecol="DSS.time",subset.mat = NULL, show.conf.int = F,
           contcol = g,title= paste("TCGA ",g),percentile.colors=c("blue","grey","red",'black'),
           time.mult=365,use.cont = F,ylab="Desease Specific Survival",max.x = 15,print.significant.only = F)
  
  
  
  
}
dev.off()




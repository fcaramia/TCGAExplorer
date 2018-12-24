rm(list = ls())
library(data.table)
library(openxlsx)
library(dplyr)
library(genefu)
source("~/Documents/Work/Projects/Sherene/Combined_Plotting/survival_functions.R")
source('PropensityScoresFunctions.R')


output_dir="~/Documents/PhD/GenderAnalysis/TCGA/Analysis/TCGAExpressionExplorerOutput/SurvivalAnalysis/"
if (!file.exists(output_dir)){
  dir.create(file.path(output_dir))
}


expression = fread("~/Documents/PhD/Data/TCGA_Xena/Expression/EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena")
datasets = c("ESCA","HNSC","LUSC","BLCA","LIHC","STAD","LGG","COAD","PAAD","READ","SKCM",'LUAD')
clinical = read.xlsx("~/Documents/PhD/Data/TCGA_CLINICAL/mmc1.xlsx", sheet = 1)
clinical = filter(clinical, type%in%datasets)
clinical$PatientID = gsub("TCGA\\-([[:alnum:]]{2})\\-([[:alnum:]]{4}).*","\\2",clinical$bcr_patient_barcode)
##P53 Interactors
p53.int = read.csv("~/Documents/PhD/GenderAnalysis/TP53_interactions/Candidates.2018.12.20.full.csv")
p53.low = as.vector(p53.int[which(p53.int$experimentally_determined_interaction>=0.3),'GeneSymbol'])
p53.low = as.vector(p53.int[,'GeneSymbol'])

####Gene annotation for EntrezGeneID
gene.annot = read.csv("~/Documents/PhD/GenderAnalysis/TCGA/Analysis/annotation.all.genes.csv",as.is = T)
gene.annot = gene.annot[,-1]
gene.annot = gene.annot[gene.annot$Chrm!='',]
gene.annot = gene.annot[order(gene.annot$Chrm,gene.annot$Entrez.Gene.ID),]

p53.int$EntrezID = gene.annot[match(p53.int$GeneSymbol,gene.annot$Approved.Symbol),'Entrez.Gene.ID']

#####Filter p53 genes
expression <- subset(expression, select=which(!duplicated(names(expression)))) 
expression = filter(expression,sample%in%p53.int$GeneSymbol)
rownames(expression) = expression$sample
expression$sample = NULL
expression = as.data.frame(t(expression))
expression$Sample = rownames(expression)
expression$SampleType = gsub("TCGA\\-([[:alnum:]]{2})\\-([[:alnum:]]{4})\\-([[:alnum:]]{2})","\\3",rownames(expression))
expression$PatientID = gsub("TCGA\\-([[:alnum:]]{2})\\-([[:alnum:]]{4}).*","\\2",rownames(expression))
expression = filter(expression,PatientID%in%clinical$PatientID,SampleType!='11')
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

mat$p53.status = ifelse(mat$PatientID%in%p53.mutants,'Mt-p53','Wt-p53')



# ####Discretize p53 string mutations by patient
DNA.muts = fread("~/Documents/PhD/GenderAnalysis/TCGA/Analysis/full.reduced.all.raw.TCGA.curated.mutations.csv")
DNA.muts.p53.string = filter(DNA.muts,HUGO_SYMBOL%in%as.vector(sig$probe ))
as.data.frame(DNA.muts.p53.string %>% group_by(PATIENT_ID) %>% dplyr::summarise(N.P = n())) -> desc.p53.string.muts.numbers
mat$p53.string.muts.desc = desc.p53.string.muts.numbers[match(mat$PatientID,desc.p53.string.muts.numbers$PATIENT_ID),'N.P']
mat[which(is.na(mat$p53.string.muts.desc)),'p53.string.muts.desc'] = 0

mat$p53.string.muts.desc.d = ifelse(mat$p53.string.muts.desc==0,'Wt',ifelse(mat$p53.string.muts.desc>6,NA,ifelse(mat$p53.string.muts.desc>2,'3 to 6','1 to 2')))

pdf(paste(output_dir,"/TCGA_p53_set_mutations_survival_os.pdf",sep=""))

mat$comb = paste(mat$p53.string.muts.desc.d,mat$p53.status,sep = '.')
mat2 = mat[!is.na(mat$p53.string.muts.desc.d)&mat$p53.string.muts.desc.d!='Wt',]

survplot(mat2,quantcol="p53.string.muts.desc.d",survcol="OS",timecol="OS.time",subset.mat = NULL, show.conf.int = F,
         contcol = 'p53.string.muts.desc.d',title= paste("TCGA p53 string muts"),percentile.colors=c("blue","grey","red",'black'),
         time.mult=365,use.cont = F,ylab="Overall Survival",max.x = 5)


survplot(mat2,quantcol="comb",survcol="OS",timecol="OS.time",subset.mat = NULL, show.conf.int = F,
         contcol = 'comb',title= paste("TCGA Gender"),percentile.colors=c("blue","grey","red",'black','yellow','green','orange','purple'),
         time.mult=365,use.cont = F,ylab="Overall Survival",max.x = 5)


dev.off()



####Do plots


genes = c(as.vector(sig$probe),'p53.string.sig')
#genes = c('p53.string.sig')

###Male and female 


pdf(paste(output_dir,"/TCGA_gender_survival_os.pdf",sep=""))
survplot(mat,quantcol="gender",survcol="OS",timecol="OS.time",subset.mat = NULL, show.conf.int = F,
         contcol = 'gender',title= paste("TCGA Gender"),percentile.colors=c("blue","grey","red",'black'),
         time.mult=365,use.cont = F,ylab="Overall Survival",max.x = 15)
mat$comb = paste(mat$gender,mat$p53.status,sep = '.')
survplot(mat,quantcol="comb",survcol="OS",timecol="OS.time",subset.mat = NULL, show.conf.int = F,
         contcol = 'comb',title= paste("TCGA Gender"),percentile.colors=c("blue","grey","red",'black'),
         time.mult=365,use.cont = F,ylab="Overall Survival",max.x = 15)


dev.off()

# for (g in genes){
# 
#   gene.d = paste(g,'d',sep = '.')
#   mat[,gene.d] = ifelse(mat[,g]>=median(mat[,g],na.rm = T),paste('hi',g,sep = '-'),paste('lo',g,sep = '-'))
#   mat$comb = paste(mat$gender,mat[,gene.d],sep = '.')
# 
#   pdf(paste(output_dir,"/TCGA_",g,"_gender_survival_os.pdf",sep=""))
#   mat2 = mat[which(!is.na(mat[,g])),]
#   survplot(mat2,quantcol="comb",survcol="OS",timecol="OS.time",subset.mat = NULL, show.conf.int = F,
#            contcol = g,title= paste("TCGA ",g),percentile.colors=c("blue","grey","red",'black'),
#            time.mult=365,use.cont = F,ylab="Overall Survival",max.x = 15)
# 
#   dev.off()
# 
# 
# }

#####p53 Wt vs p53 Mt
###For Sue

genes = c(as.vector(sig$probe),'p53.string.sig')
#genes = c('p53.string.sig')

pdf(paste(output_dir,"/TCGA_Mtp53_status_survival_25-75prct.pdf",sep=""))

for (g in genes){

  gene.d = paste(g,'d',sep = '.')
  q.sig = quantile(mat[,g],probs = c(.25,.75),na.rm = T)
  mat[,gene.d] = ifelse(mat[,g]>=q.sig[2],paste('hi',g,sep = '-'),
                           ifelse(mat[,g]<=q.sig[1],paste('lo',g,sep = '-'),NA))
  #mat[,gene.d] = ifelse(mat[,g]>=median(mat[,g],na.rm = T),paste('hi',g,sep = '-'),paste('lo',g,sep = '-'))

  mat$comb = paste(mat$p53.status,mat[,gene.d],sep = '.')


  mat2 = mat[which(!is.na(mat[,g])&!is.na(mat[,gene.d])&mat[,"p53.status"]=='Mt-p53'),]
  #mat2=mat[which(!is.na(mat[,g])&!is.na(mat[,gene.d])),]
  survplot(mat2,quantcol="comb",survcol="DSS",timecol="DSS.time",subset.mat = NULL, show.conf.int = F,
           contcol = g,title= paste("TCGA ",g),percentile.colors=c("blue","grey","red",'black'),
           time.mult=365,use.cont = F,ylab="Desease Specific Survival",max.x = 5,print.significant.only = T)




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




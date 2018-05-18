rm(list = ls())
source("ExpressionPlotsFunctions.R")
library(dplyr)
#TCGA Expression plots
#Datasets and directories and variables
datasets = c("ESCA","HNSC","LUSC","BLCA","LIHC","STAD","LGG","COAD","PAAD","READ","SKCM")
#datasets = c("READ")
output.dir = "~/Documents/PhD/GenderAnalysis/TCGA/Analysis/TCGAExpressionExplorerOutput/"

###Mutations####
mutations = read.csv("~/Documents/PhD/GenderAnalysis/TCGA/Analysis/full.reduced.all.raw.TCGA.curated.mutations.csv", as.is = T)
mutations$PATIENT_ID = gsub("TCGA\\-([[:alnum:]]{2})\\-([[:alnum:]]{4}).*","\\2", mutations$SAMPLE_ID)
###############

##############
###Read DNA VAF#####
mutations %>% filter(CANCER_TYPE%in%datasets) -> mutations
mutations %>% group_by(GENDER,CANCER_TYPE,PATIENT_ID) %>% dplyr::summarise() %>% 
  group_by(GENDER,CANCER_TYPE) %>% dplyr::summarise(N.P = n()) -> patient.numbers.by.cancer


mutations %>% group_by(GENDER,PATIENT_ID) %>% dplyr::summarise() %>% 
  group_by(GENDER) %>% dplyr::summarise(N.P = n()) -> patient.numbers


#############


filter.out = c("Silent",'Intron','IGR',"In_Frame_Ins" ,"In_Frame_Del", "lincRNA" )

filter(mutations,!(VARIANT_CLASSIFICATION%in%filter.out)) -> mutations

###Clinical Data####
demo = read.csv("~/Documents/PhD/GenderAnalysis/TCGA/Analysis/all.TCGA.curated.clinical.csv", as.is = T)
demo$PATIENT_ID = toupper(demo$PATIENT_ID)
###############

#READ TP53 Interactors

p53.int = read.csv("~/Documents/PhD/GenderAnalysis/TP53_interactions/Candidates.2018.03.01.csv")

p53.high = as.vector(p53.int[which(p53.int$score>=0.6),'ID'])
p53.med = as.vector(p53.int[which(p53.int$score>=0.4),'ID'])
p53.low = as.vector(p53.int$ID)


###Contrast VARS###
comp.clinical.vars = c("GENDER")

##Summarise mutations
mutations %>% 
  left_join(demo[c('GENDER','PATIENT_ID')]) %>%
  group_by(PATIENT_ID,GENDER,CANCER_TYPE) %>%
  dplyr::summarise(nmuts=log2(n())) %>% 
  filter(!is.na(GENDER)) -> mutations.per.patient

mutations %>% filter(HUGO_SYMBOL%in%p53.high) %>%
  left_join(demo[c('GENDER','PATIENT_ID')]) %>%
  group_by(PATIENT_ID,GENDER,CANCER_TYPE) %>%
  dplyr::summarise(nmuts=log2(n())) %>% 
  filter(!is.na(GENDER)) -> mutations.per.patient.p53.high
mutations %>% filter(HUGO_SYMBOL%in%p53.med) %>%
  left_join(demo[c('GENDER','PATIENT_ID')]) %>%
  group_by(PATIENT_ID,GENDER,CANCER_TYPE) %>%
  dplyr::summarise(nmuts=log2(n())) %>% 
  filter(!is.na(GENDER)) -> mutations.per.patient.p53.med
mutations %>% filter(HUGO_SYMBOL%in%p53.low) %>%
  left_join(demo[c('GENDER','PATIENT_ID')]) %>%
  group_by(PATIENT_ID,GENDER,CANCER_TYPE) %>%
  dplyr::summarise(nmuts=log2(n())) %>% 
  filter(!is.na(GENDER)) -> mutations.per.patient.p53.low




mutations %>% 
  filter(HUGO_SYMBOL=='TP53') %>%
  select(PATIENT_ID) -> mutants

mutations %>% mutate(TP53_STATUS = ifelse(PATIENT_ID%in%mutants$PATIENT_ID,'mt','wt')) -> mutations

mutations.per.patient$GENDER = as.factor(mutations.per.patient$GENDER)

##All Chromosomes

pvalue <- try(summary(lm(mutations.per.patient$nmuts~mutations.per.patient$GENDER))$coef[2,4])

pvalue.all = format(round(pvalue, 5), nsmall = 5)
nmales = dim(mutations.per.patient[which(mutations.per.patient$GENDER=='male'),])[1]
nfemales = dim(mutations.per.patient[which(mutations.per.patient$GENDER=='female'),])[1]   

#The X chrom

mutations %>% 
  left_join(demo[c('GENDER','PATIENT_ID')]) %>%
  filter(CHROMOSOME=='X') %>%
  group_by(PATIENT_ID,GENDER,CANCER_TYPE, CHROMOSOME) %>%
  dplyr::summarise(nmuts=log2(n())) %>% 
  filter(!is.na(GENDER)) -> mutations.per.patient.x

mutations.per.patient.x$GENDER = as.factor(mutations.per.patient.x$GENDER)

pvalue <- try(summary(lm(mutations.per.patient.x$nmuts~mutations.per.patient.x$GENDER))$coef[2,4])

pvalue.x = format(round(pvalue, 5), nsmall = 5)
nmales.x = dim(mutations.per.patient[which(mutations.per.patient.x$GENDER=='male'),])[1]
nfemales.x = dim(mutations.per.patient[which(mutations.per.patient.x$GENDER=='female'),])[1] 
###Do the plots
dir.create(paste(output.dir,'AllDataSets',"/MutationPlots",sep = ""),showWarnings = F)

tiff(paste(output.dir,'AllDataSets',"/MutationPlots/",'ExomeMutationPlotAllChrom.tiff',sep = ""),width = 400, height = 400,res = 125)

DoMutNumPlotGender(data = mutations.per.patient, data.col = 'nmuts' ,
                   fill.col = 'GENDER', x.lab = 'GENDER', 
                    y.lab = 'No. of exome mutations per patient (log2)',
                   title.txt = paste('Disparity-set: Autosomes &\n X chromosome'), nmales = 2719, nfemales = 1458)
dev.off()





# DoMutNumPlotGenderFacet(data = mutations.per.patient, data.col = 'nmuts' , facet.col = 'CANCER_TYPE',
#                    fill.col = 'GENDER', x.lab = 'Gender', 
#                    y.lab = 'log2 Mutation Numbers',
#                    title.txt = paste('Mutation Numbers, All Cancers'))

tiff(paste(output.dir,'AllDataSets',"/MutationPlots/",'ExomeMutationPlotXChrom.tiff',sep = ""),width = 400, height = 400,res=125)

DoMutNumPlotGender(data = mutations.per.patient.x, data.col = 'nmuts' ,
                   fill.col = 'GENDER', x.lab = 'GENDER', 
                   y.lab = 'No. of exome mutations per patient (log2)',
                   title.txt = paste('Disparity-set:\nX chromosome'), nmales = 2719, nfemales = 1458)



dev.off()

tiff(paste(output.dir,'AllDataSets',"/MutationPlots/",'ExomeMutationPlotXChromAllCancersXChrom.p53.low.tiff',sep = ""),width = 600, height = 400,res=125)

 DoMutNumPlotGenderFacet(data = mutations.per.patient.p53.low, data.col = 'nmuts' , facet.col = 'CANCER_TYPE',
                         fill.col = 'GENDER', x.lab = 'Gender', 
                         y.lab = 'log2 Mutation Numbers',
                         title.txt = paste('Mutation Numbers, All Cancers, X Chromosome'))

 dev.off()
 
#MUTANTS

mutations %>% 
  filter(TP53_STATUS == 'mt') %>%
  left_join(demo[c('GENDER','PATIENT_ID')]) %>%
  group_by(PATIENT_ID,GENDER,CANCER_TYPE) %>%
  dplyr::summarise(nmuts=log2(n())) %>% 
  filter(!is.na(GENDER)) -> mutations.per.patient

jpeg(paste(output.dir,'AllDataSets',"/MutationPlots/",'ExomeMutationPlotAllChromp53Mutants.jpeg',sep = ""),width = 800, height = 1200,res=300)

DoMutNumPlotGender(data = mutations.per.patient, data.col = 'nmuts' ,
                   fill.col = 'GENDER', x.lab = 'p53 Mutants',
                   y.lab = 'No. of exome mutations per patient (log2)',
                   title.txt = paste('Disparity-set: Autosomes &\n X chromosome\np53 Mutants'))
dev.off()


mutations %>% 
  filter(TP53_STATUS == 'wt') %>%
  left_join(demo[c('GENDER','PATIENT_ID')]) %>%
  group_by(PATIENT_ID,GENDER,CANCER_TYPE) %>%
  dplyr::summarise(nmuts=log2(n())) %>% 
  filter(!is.na(GENDER)) -> mutations.per.patient

jpeg(paste(output.dir,'AllDataSets',"/MutationPlots/",'ExomeMutationPlotAllChromp53Wt.jpeg',sep = ""),width = 800, height = 1200,res=300)

DoMutNumPlotGender(data = mutations.per.patient, data.col = 'nmuts' ,
                   fill.col = 'GENDER', x.lab = 'p53 Wt',
                   y.lab = 'No. of exome mutations per patient (log2)',
                   title.txt = paste('Disparity-set: Autosomes &\n X chromosome\np53 Wt'))
dev.off()

###X CHROM



mutations %>% 
  filter(TP53_STATUS == 'mt' & CHROMOSOME=='X') %>%
  left_join(demo[c('GENDER','PATIENT_ID')]) %>%
  group_by(PATIENT_ID,GENDER,CANCER_TYPE) %>%
  dplyr::summarise(nmuts=log2(n())) %>% 
  filter(!is.na(GENDER)) -> mutations.per.patient

jpeg(paste(output.dir,'AllDataSets',"/MutationPlots/",'ExomeMutationPlotXChromp53Mutants.jpeg',sep = ""),width = 800, height = 1200,res=300)

DoMutNumPlotGender(data = mutations.per.patient, data.col = 'nmuts' ,
                   fill.col = 'GENDER', x.lab = 'p53 Mutants',
                   y.lab = 'No. of exome mutations per patient (log2)',
                   title.txt = paste('Disparity-set: X chromosome\np53 Mutants'))
dev.off()


mutations %>% 
  filter(TP53_STATUS == 'wt' & CHROMOSOME=='X') %>%
  left_join(demo[c('GENDER','PATIENT_ID')]) %>%
  group_by(PATIENT_ID,GENDER,CANCER_TYPE) %>%
  dplyr::summarise(nmuts=log2(n())) %>% 
  filter(!is.na(GENDER)) -> mutations.per.patient

jpeg(paste(output.dir,'AllDataSets',"/MutationPlots/",'ExomeMutationPlotXChromp53Wt.jpeg',sep = ""),width = 800, height = 1200,res=300)

DoMutNumPlotGender(data = mutations.per.patient, data.col = 'nmuts' ,
                   fill.col = 'GENDER', x.lab = 'p53 Wt',
                   y.lab = 'No. of exome mutations per patient (log2)',
                   title.txt = paste('Disparity-set: X chromosome\np53 Wt'))
dev.off()

###COMPARE MUTANTS to WT

mutations %>% 
  left_join(demo[c('GENDER','PATIENT_ID')]) %>%
  group_by(PATIENT_ID,TP53_STATUS,CANCER_TYPE) %>%
  dplyr::summarise(nmuts=log2(n())) -> mutations.per.patient

jpeg(paste(output.dir,'AllDataSets',"/MutationPlots/",'ExomeMutationPlotAllChromMtVSWt.jpeg',sep = ""),width = 800, height = 1200,res=300)

DoMutNumPlotGender(data = mutations.per.patient, data.col = 'nmuts' ,
                   fill.col = 'TP53_STATUS', x.lab = 'p53 status',
                   y.lab = 'No. of exome mutations per patient (log2)',
                   title.txt = paste('Disparity-set: All chromosomes'))
dev.off()

mutations %>% 
  filter(CHROMOSOME=='X') %>%
  left_join(demo[c('GENDER','PATIENT_ID')]) %>%
  group_by(PATIENT_ID,TP53_STATUS,CANCER_TYPE) %>%
  dplyr::summarise(nmuts=log2(n())) -> mutations.per.patient

jpeg(paste(output.dir,'AllDataSets',"/MutationPlots/",'ExomeMutationPlotXChromMtVSWt.jpeg',sep = ""),width = 800, height = 1200,res=300)

DoMutNumPlotGender(data = mutations.per.patient, data.col = 'nmuts' ,
                   fill.col = 'TP53_STATUS', x.lab = 'p53 status',
                   y.lab = 'No. of exome mutations per patient (log2)',
                   title.txt = paste('Disparity-set: X chromosome'))
dev.off()



muts.bkup = mutations
###Start Main Loop###
for (i in datasets)
{
  print(i)
  ###Create directories####
  dir.create(paste(output.dir,i,sep = ""),showWarnings = F)
  dir.create(paste(output.dir,i,"/MutationPlots/",sep = ""),showWarnings = F)
  #########################
  
  muts.bkup %>%
    filter(CANCER_TYPE == i) -> mutations
  
  ##Summarise mutations
  mutations %>% 
    left_join(demo[c('GENDER','PATIENT_ID')]) %>%
    group_by(PATIENT_ID,GENDER,CANCER_TYPE) %>%
    dplyr::summarise(nmuts=log2(n())) %>% 
    filter(!is.na(GENDER)) -> mutations.per.patient
  
  mutations.per.patient$GENDER = as.factor(mutations.per.patient$GENDER)
  
  pvalue <- try(summary(lm(mutations.per.patient$nmuts~mutations.per.patient$GENDER))$coef[2,4])
  
  pvalue.all = format(round(pvalue, 5), nsmall = 5)
  nmales = dim(mutations.per.patient[which(mutations.per.patient$GENDER=='male'),])[1]
  nfemales = dim(mutations.per.patient[which(mutations.per.patient$GENDER=='female'),])[1]   
  
  #The X chrom
  
  mutations %>% 
    filter(CHROMOSOME=='X') %>%
    left_join(demo[c('GENDER','PATIENT_ID')]) %>%
    group_by(PATIENT_ID,GENDER,CANCER_TYPE) %>%
    dplyr::summarise(nmuts=log2(n())) %>% 
    filter(!is.na(GENDER)) ->mutations.per.patient.x
  
  
  mutations.per.patient.x$GENDER = as.factor(mutations.per.patient.x$GENDER)
  
  pvalue <- try(summary(lm(mutations.per.patient.x$nmuts~mutations.per.patient.x$GENDER))$coef[2,4])
  
  pvalue.x = format(round(pvalue, 5), nsmall = 5)
  nmales.x = dim(mutations.per.patient[which(mutations.per.patient.x$GENDER=='male'),])[1]
  nfemales.x = dim(mutations.per.patient[which(mutations.per.patient.x$GENDER=='female'),])[1] 
  #####PLOTS
  
  jpeg(paste(output.dir,i,"/MutationPlots/",'ExomeMutationPlotAllChrom.jpeg',sep = ""),width = 800, height = 1200,res=300)
  

    
  DoMutNumPlotGender(data = mutations.per.patient, data.col = 'nmuts' ,
                     fill.col = 'GENDER', x.lab = paste('p-val = ',pvalue.all), 
                     y.lab = 'No. of exome mutations per patient (log2)',
                     title.txt = paste(i,': Autosomes & \nX chromosome'))
  dev.off()
  
  jpeg(paste(output.dir,i,"/MutationPlots/",'ExomeMutationPlotXChrom.jpeg',sep = ""),width = 800, height = 1200,res=300)
  
  DoMutNumPlotGender(data = mutations.per.patient.x, data.col = 'nmuts' ,
                     fill.col = 'GENDER', x.lab = paste('p-val = ',pvalue.x), 
                     y.lab = 'No. of exome mutations per patient (log2)',
                     title.txt = paste(i,': \nX chromosome'))
  
  
  dev.off()
  
  
  #MUTANTS
  
  mutations %>% 
    filter(TP53_STATUS == 'mt') %>%
    left_join(demo[c('GENDER','PATIENT_ID')]) %>%
    group_by(PATIENT_ID,GENDER,CANCER_TYPE) %>%
    dplyr::summarise(nmuts=log2(n())) %>% 
    filter(!is.na(GENDER)) -> mutations.per.patient
  
  jpeg(paste(output.dir,i,"/MutationPlots/",'ExomeMutationPlotAllChromp53Mutants.jpeg',sep = ""),width = 800, height = 1200,res=300)
  
  DoMutNumPlotGender(data = mutations.per.patient, data.col = 'nmuts' ,
                     fill.col = 'GENDER', x.lab = 'p53 Mutants',
                     y.lab = 'No. of exome mutations per patient (log2)',
                     title.txt = paste('Disparity-set: Autosomes &\n X chromosome\np53 Mutants'))
  dev.off()
  
  
  mutations %>% 
    filter(TP53_STATUS == 'wt') %>%
    left_join(demo[c('GENDER','PATIENT_ID')]) %>%
    group_by(PATIENT_ID,GENDER,CANCER_TYPE) %>%
    dplyr::summarise(nmuts=log2(n())) %>% 
    filter(!is.na(GENDER)) -> mutations.per.patient
  
  jpeg(paste(output.dir,i,"/MutationPlots/",'ExomeMutationPlotAllChromp53Wt.jpeg',sep = ""),width = 800, height = 1200,res=300)
  
  DoMutNumPlotGender(data = mutations.per.patient, data.col = 'nmuts' ,
                     fill.col = 'GENDER', x.lab = 'p53 Wt',
                     y.lab = 'No. of exome mutations per patient (log2)',
                     title.txt = paste('Disparity-set: Autosomes &\n X chromosome\np53 Wt'))
  dev.off()
  
  ###X CHROM
  
  
  
  mutations %>% 
    filter(TP53_STATUS == 'mt' & CHROMOSOME=='X') %>%
    left_join(demo[c('GENDER','PATIENT_ID')]) %>%
    group_by(PATIENT_ID,GENDER,CANCER_TYPE) %>%
    dplyr::summarise(nmuts=log2(n())) %>% 
    filter(!is.na(GENDER)) -> mutations.per.patient
  
  jpeg(paste(output.dir,i,"/MutationPlots/",'ExomeMutationPlotXChromp53Mutants.jpeg',sep = ""),width = 800, height = 1200,res=300)
  
  DoMutNumPlotGender(data = mutations.per.patient, data.col = 'nmuts' ,
                     fill.col = 'GENDER', x.lab = 'p53 Mutants',
                     y.lab = 'No. of exome mutations per patient (log2)',
                     title.txt = paste('Disparity-set: X chromosome\np53 Mutants'))
  dev.off()
  
  
  mutations %>% 
    filter(TP53_STATUS == 'wt' & CHROMOSOME=='X') %>%
    left_join(demo[c('GENDER','PATIENT_ID')]) %>%
    group_by(PATIENT_ID,GENDER,CANCER_TYPE) %>%
    dplyr::summarise(nmuts=log2(n())) %>% 
    filter(!is.na(GENDER)) -> mutations.per.patient
  
  jpeg(paste(output.dir,i,"/MutationPlots/",'ExomeMutationPlotXChromp53Wt.jpeg',sep = ""),width = 800, height = 1200,res=300)
  
  DoMutNumPlotGender(data = mutations.per.patient, data.col = 'nmuts' ,
                     fill.col = 'GENDER', x.lab = 'p53 Wt',
                     y.lab = 'No. of exome mutations per patient (log2)',
                     title.txt = paste('Disparity-set: X chromosome\np53 Wt'))
  dev.off()
  
  ###COMPARE MUTANTS to WT
  
  mutations %>% 
    left_join(demo[c('GENDER','PATIENT_ID')]) %>%
    group_by(PATIENT_ID,TP53_STATUS,CANCER_TYPE) %>%
    dplyr::summarise(nmuts=log2(n())) -> mutations.per.patient
  
  jpeg(paste(output.dir,i,"/MutationPlots/",'ExomeMutationPlotAllChromMtVSWt.jpeg',sep = ""),width = 800, height = 1200,res=300)
  
  DoMutNumPlotGender(data = mutations.per.patient, data.col = 'nmuts' ,
                     fill.col = 'TP53_STATUS', x.lab = 'p53 status',
                     y.lab = 'No. of exome mutations per patient (log2)',
                     title.txt = paste('Disparity-set: All chromosomes'))
  dev.off()
  
  mutations %>% 
    filter(CHROMOSOME=='X') %>%
    left_join(demo[c('GENDER','PATIENT_ID')]) %>%
    group_by(PATIENT_ID,TP53_STATUS,CANCER_TYPE) %>%
    dplyr::summarise(nmuts=log2(n())) -> mutations.per.patient
  
  jpeg(paste(output.dir,i,"/MutationPlots/",'ExomeMutationPlotXChromMtVSWt.jpeg',sep = ""),width = 800, height = 1200,res=300)
  
  DoMutNumPlotGender(data = mutations.per.patient, data.col = 'nmuts' ,
                     fill.col = 'TP53_STATUS', x.lab = 'p53 status',
                     y.lab = 'No. of exome mutations per patient (log2)',
                     title.txt = paste('Disparity-set: X chromosome'))
  dev.off()
  
  
  
}

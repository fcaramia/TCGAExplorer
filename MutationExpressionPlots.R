rm(list = ls())
library(dplyr)
library(ggplot2)
library(gtools)
source('ExpressionPlotsFunctions.R')
#TCGA Expression plots
#Datasets and directories and variables
datasets = c("ESCA","HNSC","LUSC","BLCA","LIHC","STAD","LGG","COAD","PAAD","READ","SKCM")
output.dir = "~/Documents/PhD/GenderAnalysis/TCGA/Analysis/TCGAExpressionExplorerOutput/"

##RMAF DATA
RMAF = read.csv("~/Documents/PhD/GenderAnalysis/TCGA/Analysis/all.TCGA.ExMut.V2.csv", as.is=T)
RMAF$PATIENT_ID = gsub("TCGA\\-([[:alnum:]]{2})\\-([[:alnum:]]{4}).*","\\2", RMAF$SAMPLE)

signature.file =  "~/Documents/PhD/GenderAnalysis/TCGA/Analysis/TP53.lists.2.txt"
###Read Signatures####
sig.lists = read.delim(signature.file, sep = " ")

####Read mutations
mutations = read.csv("~/Documents/PhD/GenderAnalysis/TCGA/Analysis/reduced.all.TCGA.curated.mutations.csv", as.is = T)
mutations$PATIENT_ID = gsub("TCGA\\-([[:alnum:]]{2})\\-([[:alnum:]]{4}).*","\\2", mutations$SAMPLE_ID)


mutations%>% filter(HUGO_SYMBOL=='TP53')%>% select(PATIENT_ID) -> p53.mutants 

#Clean low reads
RMAF$SUM = RMAF$REF_COUNT + RMAF$ALT_COUNT
RMAF %>% filter(SUM>10&VARIANT_CLASSIFICATION!='Silent') -> RMAF
RMAF$RMAF = RMAF$ALT_COUNT/RMAF$SUM
RMAF$LOGIT.RMAF = logit(x = RMAF$RMAF,min = 0-0.1 , max = 1+0.1)
RMAF %>% mutate(P53.STATUS = ifelse(PATIENT_ID%in%p53.mutants$PATIENT_ID,'Mt','Wt')) -> RMAF


##Do something about noise?######

###Clinical Data####
demo = read.csv("~/Documents/PhD/GenderAnalysis/TCGA/Analysis/all.TCGA.curated.clinical.csv", as.is = T)
demo$PATIENT_ID = toupper(demo$PATIENT_ID)
###############




pvalue <- try(summary(lm(RMAF$LOGIT.RMAF~RMAF$GENDER))$coef[2,4])

pvalue.all = ifelse(pvalue>0.0001,format(pvalue, scientific=T, digits=2),0)

RMAF%>% filter(GENDER=='female') %>% dim() -> nfemales
RMAF%>% filter(GENDER=='male') %>% dim() -> nmales

RMAF%>% filter(GENDER=='female') %>% 
  group_by(PATIENT_ID) %>% 
  dplyr::summarise(NFEM = n()) %>% 
  dim() -> females
RMAF %>% filter(GENDER=='male') %>% 
  group_by(PATIENT_ID) %>% 
  dplyr::summarise(NFEM = n()) %>% 
  dim() -> males


###Do the plots
dir.create(paste(output.dir,'AllDataSets',"/RMAFPlots",sep = ""),showWarnings = F)

png(paste(output.dir,'AllDataSets',"/RMAFPlots/",'X-CHROM.RMAFPlots.V2.png',sep = ""),width = 1200, height = 800,res=170)


DoDensityPlotGenderRMAF(dat = RMAF, dat.col = 'RMAF' ,
                   fill.col = 'GENDER', 
                   title.txt = paste('RMAF numbers x-chromosomes, All Cancers','\np-value = ',0,
                                     '\nExome mutations with RNA content in Males = ',nmales[1],
                                     '\nExome mutations with RNA content in Females = ',nfemales[1],
                                     '\nMale Patients analysed = ', males[1],
                                     '\nFemale Patients analysed = ', females[1]))
dev.off()

png(paste(output.dir,'AllDataSets',"/RMAFPlots/",'X-CHROM.RMAFPlots.P53.V2.png',sep = ""),width = 1200, height = 800,res=170)

DoDensityPlotFacetGenderRMAF(dat = RMAF, dat.col = 'RMAF', facet.col = 'P53.STATUS'
                             ,fill.col = 'GENDER',
                             title.txt = paste("X-CHROM",
                                               'RMAF densities by Gender and p-53 status'))  

dev.off()


png(paste(output.dir,'AllDataSets',"/RMAFPlots/",'X-CHROM.by.Cancer.RMAFPlots.V2.png',sep = ""),width = 1200, height = 800, res =120)
RMAF %>% 
  group_by(CANCER_TYPE) %>% 
  dplyr::mutate(NMUTS = n()) -> 
  RMAF
RMAF$LABEL = paste(RMAF$CANCER_TYPE,'MUTS:',RMAF$NMUTS)
DoDensityPlotFacetGenderRMAF(dat = RMAF, dat.col = 'RMAF', facet.col = 'LABEL'
                              ,fill.col = 'GENDER',title.txt = 'RMAF densities by Cancer')
dev.off()




for (s in rownames(sig.lists)){
###############by signature
  #Add signature
  genes = unlist(strsplit(as.character(sig.lists[s,'GENES']),',')) 
  
  RMAF %>%
    left_join(demo[c('GENDER','PATIENT_ID')]) %>%
    #left_join(weights.all) %>%
    filter(HUGO_SYMBOL%in%genes & !is.na(GENDER)) -> RMAF.sig
  
  
  
  pvalue <- try(summary(lm(RMAF.sig$LOGIT.RMAF~RMAF.sig$GENDER))$coef[2,4])
  
  pvalue.all = ifelse(pvalue>0.0001,format(pvalue, scientific=T, digits=2),0)
  
  RMAF.sig %>% filter(GENDER=='female') %>% dim() -> nfemales
  RMAF.sig %>% filter(GENDER=='male') %>% dim() -> nmales
  
  RMAF.sig %>% filter(GENDER=='female') %>% 
    group_by(PATIENT_ID) %>% 
    dplyr::summarise(NFEM = n()) %>% 
    dim() -> females
  RMAF.sig %>% filter(GENDER=='male') %>% 
    group_by(PATIENT_ID) %>% 
    dplyr::summarise(NFEM = n()) %>% 
    dim() -> males
  
  RMAF.sig %>% 
    group_by(CANCER_TYPE) %>% 
    dplyr::mutate(NMUTS = n()) -> 
    RMAF.sig
  RMAF.sig$LABEL = paste(RMAF.sig$CANCER_TYPE,'MUTS:',RMAF.sig$NMUTS)
  
  ###Do the plots
  
  png(paste(output.dir,'AllDataSets',"/RMAFPlots/",as.character(sig.lists[s,'SIGNATURE']),'_RMAFPlot.V2.png',sep = ""),width = 1200, height = 800,res=170)
  
  DoDensityPlotGenderRMAF(dat = RMAF.sig, dat.col = 'LOGIT.RMAF' ,
                          fill.col = 'GENDER', 
                          title.txt = paste('RMAF numbers X-chromosome, All Cancers','\np-value = ',0,
                                            '\nExome mutations with RNA content in Males = ',nmales[1],
                                            '\nExome mutations with RNA content in Females = ',nfemales[1],
                                            '\nMale Patients analysed = ', males[1],
                                            '\nFemale Patients analysed = ', females[1]))
  dev.off()

  
  png(paste(output.dir,'AllDataSets',"/RMAFPlots/",as.character(sig.lists[s,'SIGNATURE']),'_RMAFPlotbyCancer.P53.V2.png',sep = ""),width = 1200, height = 800, res=120)
  
  
  RMAF.sig %>% filter(GENDER=='female'&P53.STATUS == 'Mt') %>% dim() -> nfemales.Mt
  RMAF.sig %>% filter(GENDER=='male'&P53.STATUS == 'Mt') %>% dim() -> nmales.Mt
  DoDensityPlotFacetGenderRMAF(dat = RMAF.sig, dat.col = 'RMAF', facet.col = 'P53.STATUS'
                               ,fill.col = 'GENDER',
                               title.txt = paste(as.character(sig.lists[s,'SIGNATURE']),
                                                 'RMAF densities by Cancer',
                                                 '\nExome mutations with RNA content in Wt Males = ',nmales[1]-nmales.Mt[1],
                                                 '\nExome mutations with RNA content in Wt Females = ',nfemales[1]-nfemales.Mt[1],
                                                 '\nExome mutations with RNA content in Mt Males = ',nmales.Mt[1],
                                                 '\nExome mutations with RNA content in Mt Females = ',nfemales.Mt[1]))  
  dev.off()
}


##########################


for (i in datasets){
  
  print(i)
  ###Create directories####
  dir.create(paste(output.dir,i,sep = ""),showWarnings = F)
  dir.create(paste(output.dir,i,"/RMAFPlots/",sep = ""),showWarnings = F)
  #########################
  
  ###Propensity Scores per Tumor Type###
  #weights.tumor = read.csv(paste(output.dir,i,'/PS/PS.csv',sep = ""))
  
  RMAF %>%
    left_join(demo[c('GENDER','PATIENT_ID')]) %>%
    #left_join(weights.tumor) %>%
    filter(CHROM=='X' & !is.na(GENDER) & CANCER_TYPE==i) -> RMAF.X
  
  RMAF.X$GENDER = as.factor(RMAF.X$GENDER)
  
  pvalue <- try(summary(lm(RMAF.X$LOGIT.RMAF~RMAF.X$GENDER))$coef[2,4])
  
  pvalue.all = ifelse(pvalue>0.0001,format(pvalue, scientific=T, digits=2),0)
  
  RMAF.X %>% filter(GENDER=='female') %>% dim() -> nfemales
  RMAF.X %>% filter(GENDER=='male') %>% dim() -> nmales
  RMAF.X %>% filter(GENDER=='female') %>% 
    group_by(PATIENT_ID) %>% 
    dplyr::summarise(NFEM = n()) %>% 
    dim() -> females
  RMAF.X %>% filter(GENDER=='male') %>% 
    group_by(PATIENT_ID) %>% 
    dplyr::summarise(NFEM = n()) %>% 
    dim() -> males
  #Do the plots
  jpeg(paste(output.dir,i,"/RMAFPlots/",'RMAFPlots.jpeg',sep = ""),width = 800, height = 800,res=150)
  DoDensityPlotGenderRMAF(dat = RMAF.X, dat.col = 'LOGIT.RMAF' ,
                        fill.col = 'GENDER', 
                        title.txt = paste('RMAF numbers X-chromosome, ',i,'\np-value = ',pvalue.all,
                                          '\nExome mutations with RNA content in Males = ',nmales[1],
                                          '\nExome mutations with RNA content in Females = ',nfemales[1],
                                          '\nMale Patients analysed = ', males[1],
                                          '\nFemale Patients analysed = ', females[1]))
  dev.off()
}

# #Filter X and Y
# #exMut = MutExpDB[which(MutExpDB$sum>=10&MutExpDB$CHROM!='X'&MutExpDB$CHROM!='Y'),]
# #silMut = MutExpDB[which(MutExpDB$sum<10&MutExpDB$CHROM!='X'&MutExpDB$CHROM!='Y'),]
# exMut = MutExpDB[which(MutExpDB$sum>=10&MutExpDB$CHROM=='X'),]
# silMut = MutExpDB[which(MutExpDB$sum<10&MutExpDB$CHROM=='X'),]
# #Calculate Percentages
# exMut$exmut.prct = exMut$ALT_COUNT/exMut$sum
# muts_per_patient = ddply(MutExpDB[which(!is.na(MutExpDB$sum)),], ~ PATIENT_ID,  nrow)
# sil_muts_per_patient = ddply(silMut, ~ PATIENT_ID,  nrow)
# 
# muts_per_patient$sil_muts = sil_muts_per_patient[match(muts_per_patient$PATIENT_ID,sil_muts_per_patient$PATIENT_ID),'V1']
# muts_per_patient[which(is.na(muts_per_patient$sil_muts)),'sil_muts'] = 0
# muts_per_patient$Gender = demo[match(muts_per_patient$PATIENT_ID,demo$PATIENT_ID),"GENDER"]
# muts_per_patient$silmut.prct = muts_per_patient$sil_muts/muts_per_patient$V1
# 
# ###Discretice Mutation expression #####
# exMut$exMut=ifelse(exMut$exmut.prct<.33, 'low', ifelse(exMut$exmut.prct<.66,"med",'high'))
# expressed_per_patient = ddply(exMut, ~ PATIENT_ID,  nrow)
# exMut_per_patient = ddply(exMut, ~ PATIENT_ID + exMut,  nrow)
# exMut_per_patient$expr.muts = expressed_per_patient[match(exMut_per_patient$PATIENT_ID,expressed_per_patient$PATIENT_ID),'V1']
# exMut_per_patient$expr.muts.d.prct = exMut_per_patient$V1/exMut_per_patient$expr.muts
# exMut_per_patient$Gender = demo[match(exMut_per_patient$PATIENT_ID,demo$PATIENT_ID),"GENDER"]
# boxplot(expr.muts.d.prct~Gender+exMut,las=1, col=c('pink',"blue"),
#         data = exMut_per_patient, 
#         main=paste("Expressed mutation per patient as %"), outline= F)
# 
# ####Silenced mutation plots######
# boxplot(silmut.prct~Gender,las=1, col=c('pink',"blue"),
#         data = muts_per_patient, 
#         main=paste("Silenced mutation per patient as %"), outline= F)
# 
# ggplot(muts_per_patient,aes(y=silmut.prct,x=Gender, fill = Gender)) +
#   #stat_boxplot(geom ='errorbar',width = 0.2, colour=c("hotpink3", "blue")) + 
#   geom_violin( ) + 
#   labs(x='Gender',y='Mutation Expression') + ggtitle(label = 'Silenced mutation per patient as %') +
#   scale_fill_manual(values=c( "lightpink",'skyblue')) + theme_minimal() 
# 
# 
# muts_per_patient$Gender.factor = as.factor(muts_per_patient$Gender)
# t.test(silmut.prct~Gender.factor,data = muts_per_patient)
# 
# ###Discretice Extreme Mutation expression #####
# exMut$exMut=ifelse(exMut$exmut.prct<.10, 'very.low', ifelse(exMut$exmut.prct<.90,"med",'very.high'))
# expressed_per_patient = ddply(exMut, ~ PATIENT_ID,  nrow)
# exMut_per_patient = ddply(exMut, ~ PATIENT_ID + exMut,  nrow)
# exMut_per_patient$expr.muts = expressed_per_patient[match(exMut_per_patient$PATIENT_ID,expressed_per_patient$PATIENT_ID),'V1']
# exMut_per_patient$expr.muts.d.prct = exMut_per_patient$V1/exMut_per_patient$expr.muts
# exMut_per_patient$Gender = demo[match(exMut_per_patient$PATIENT_ID,demo$PATIENT_ID),"GENDER"]
# boxplot(expr.muts.d.prct~Gender+exMut,las=1, col=c('pink',"blue"),
#         data = exMut_per_patient, 
#         main=paste("Expressed mutation per patient as %"), outline= F)
# 
# #####COSMIC GENES######
# COSMIC = read.csv("~/Documents/PhD/Data/COSMIC/Census_allTue Jan 17 04-15-44 2017.csv")
# boxplot(exmut.prct~Gender,las=1, col=c('pink',"blue"),
#         data = exMut[which(exMut$CHROM!='X'&exMut$HUGO_SYMBOL%in%COSMIC$Gene.Symbol),], 
#         main=paste(""), outline= F)
# 
# ggplot(exMut[which(exMut$CHROM!='X'&exMut$HUGO_SYMBOL%in%COSMIC$Gene.Symbol),],aes(y=exmut.prct,x=Gender, fill = Gender)) +
#   stat_boxplot(geom ='errorbar',width = 0.2, colour=c("hotpink3", "blue")) + 
#   geom_boxplot(outlier.colour = NA, colour=c("hotpink3", "blue")) + 
#   labs(x='Gender',y='Mutation Expression') + ggtitle(label = 'Mutation Expression COSMIC GENES 11 Cancers') + 
#   scale_fill_manual(values=c( "lightpink",'skyblue')) + theme_minimal() 
# 
# 
# exMut$Gender.factor = as.factor(exMut$Gender)
# t.test(exmut.prct~Gender.factor,data = exMut[which(exMut$CHROM!='X'&exMut$HUGO_SYMBOL%in%COSMIC$Gene.Symbol),])
# 
# 
# #All CHROM
# exMut$CHROM.G = paste(exMut$CHROM,exMut$Gender, sep = '')
# exMut = exMut[which(!is.na(exMut$Gender)),]
# ggplot(exMut,aes(y=exmut.prct,x=CHROM.G, fill=Gender)) +
#   #stat_boxplot(geom ='errorbar',width = 0.2, colour=c("hotpink3", "blue")) + 
#   geom_boxplot(outlier.colour = NA) + 
#   labs(x='Gender',y='Mutation Expression') + ggtitle(label = 'Mutation Expression All autosomes 11 Cancers') + 
#   #scale_fill_manual(values=c( "lightpink",'skyblue')) 
#    theme(axis.text.x=element_text(angle=90,hjust = 1)) 
# 
# 
# #X.chrom
# nmutation.males = nrow(exMut[exMut$Gender=='male',])
# nmutation.females =  nrow(exMut[exMut$Gender=='female',])
# 
# boxplot(exmut.prct~Gender,las=1, col=c('pink',"blue"),data = exMut[exMut$CHROM=='X',], 
#         main=paste(""), outline= F)
# 
# ggplot(exMut[which(exMut$CHROM=='X'),],aes(y=exmut.prct,x=Gender, fill = Gender)) +
#   stat_boxplot(geom ='errorbar',width = 0.2, colour=c("hotpink3", "blue")) + 
#   geom_boxplot(outlier.colour = NA, colour=c("hotpink3", "blue")) + 
#   labs(x='Gender',y='Mutation Expression') + ggtitle(label = 'Mutation Expression X Chromosome 11 Cancers') + 
#   scale_fill_manual(values=c( "lightpink",'skyblue')) + theme_minimal() 
# 
# 
# exMut$Gender.factor = as.factor(exMut$Gender)
# t.test(exmut.prct~Gender.factor,data = exMut[exMut$CHROM=='X',])
# 
# #from p53 gene list
# p53_genes = c("DKC1","FOXP3","HDAC8","LAS1L","TAF1","YY2","AMER1","BMX","CUL4B","DDX3X"
#               ,"DDX53","HUWE1","MAGEA2","MCTS1","NOX1","OGT","PSMD10","UBE2A","UTP14A",
#               "UXT","XIAP","AIFM1","G6PD","CD40L","IL2RG","SH2D1A","TLR8","AR")
# 
# exMut.gene.list = exMut[which(exMut$HUGO_SYMBOL%in%p53_genes),]
# boxplot(exmut.prct~Gender,las=1, col=c('pink',"blue"),data = exMut.gene.list)
# 
# ggplot(exMut.gene.list,aes(y=exmut.prct,x=Gender, fill = Gender)) +
#   stat_boxplot(geom ='errorbar',width = 0.2, colour=c("hotpink3", "blue")) + 
#   geom_boxplot(outlier.colour = NA, colour=c("hotpink3", "blue")) + 
#   labs(x='Gender',y='Mutation Expression') + ggtitle(label = 'Mutation Expression X Chr-p53 Linked 11 Cancers') + 
#   scale_fill_manual(values=c( "lightpink",'skyblue')) + theme_minimal() 
# 
# 
# kruskal.test(exmut.prct~Gender.factor,data = exMut.gene.list)
# 
# #####Mutation Expression per gene tables#######
# 
# aux = ddply(exMut, ~ Gender+HUGO_SYMBOL+CHROM+PATIENT_ID, summarize, max_le = max(exmut.prct))
# aux2 = ddply(aux[which(aux$max_le>0.5),], ~ Gender+HUGO_SYMBOL+CHROM,  nrow)
# write.csv(aux2,"../AllCancersStudy/plots4.0/Expressed(0.5).Mutations.by.gene.and.gender.csv", row.names = F)
# aux2 = ddply(aux[which(aux$max_le>0.9),], ~ Gender+HUGO_SYMBOL+CHROM,  nrow)
# write.csv(aux2,"../AllCancersStudy/plots4.0/Expressed(0.9).Mutations.by.gene.and.gender.csv", row.names = F)
# aux2 = ddply(aux[which(aux$max_le<=0.05),], ~ Gender+HUGO_SYMBOL+CHROM,  nrow)
# write.csv(aux2,"../AllCancersStudy/plots4.0/Expressed(less0.05).Mutations.by.gene.and.gender.csv", row.names = F)
# aux2 = ddply(aux[which(aux$max_le>0.75),], ~ Gender+HUGO_SYMBOL+CHROM,  nrow)
# write.csv(aux2,"../AllCancersStudy/plots4.0/Expressed(0.75).Mutations.by.gene.and.gender.csv", row.names = F)
# aux2 = ddply(aux, ~ Gender+HUGO_SYMBOL+CHROM,  nrow)
# aux2 = ddply(aux, ~ Gender+HUGO_SYMBOL+CHROM, summarise, med.ME = median(max_le))
# aux2 = ddply(aux, ~ Gender+HUGO_SYMBOL+CHROM, summarise, med.ME = median(max_le), n.pat=length(PATIENT_ID))
# 
# aux3 = ddply(aux2[which(aux2$n.pat>10),], ~ HUGO_SYMBOL+CHROM, summarise, diff.med.ME = diff(med.ME))
# write.csv(aux3,"~/Documents/PhD/GenderAnalysis/TCGA/Analysis/AllCancersStudy/plots4.0/Expressed.Mutations.difference.by.gender.min10.csv", row.names = F)
# write.csv(aux2,"~/Documents/PhD/GenderAnalysis/TCGA/Analysis/AllCancersStudy/plots4.0/Expressed.Mutations.medians.by.gender.csv", row.names = F)
# 
# 
# nr_males=length(unique(aux[which(aux$Gender=='male'),'PATIENT_ID']))
# nr_females=length(unique(aux[which(aux$Gender=='female'),'PATIENT_ID']))

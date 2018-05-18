rm(list = ls())
library(dplyr)
library(gtools)
source('PropensityScoresFunctions.R')
source('ExpressionPlotsFunctions.R')
Exact_Statistic <- function(a, b, c , d){
  
  n = (a+1/2)*(c+1/2)/(b+1/2)*(d+1/2)
  d = sqrt(1/(a+1/2)+1/(b+1/2)+1/(c+1/2)+1/(d+1/2))
  z = log(n)/d
  
  return(z)
}
#Datasets and directories and variables
datasets = c("ESCA","HNSC","LUSC","BLCA","LIHC","STAD","LGG","COAD","PAAD","READ","SKCM")
#datasets = c("READ")
output.dir = "~/Documents/PhD/GenderAnalysis/TCGA/Analysis/TCGAExpressionExplorerOutput/"

###Clinical Data####
clinical.data = read.csv("~/Documents/PhD/GenderAnalysis/TCGA/Analysis/all.TCGA.curated.clinical.csv", as.is = T)
clinical.data$PATIENT_ID = toupper(clinical.data$PATIENT_ID)

###############

###Propensity scores####
#Cofounding Factors definitions 
cofounding.factors.ind = c('SMOKING_STATUS','RACE','PATHOLOGIC_STAGE')
cofounding.factors.all = c('SMOKING_STATUS','RACE')
cofounding.factors.default = c('PATIENT_ID','GENDER','AGE_AT_INITIAL_PATHOLOGIC_DIAGNOSIS')
########################

#READ TP53 Interactors

p53.int = read.csv("~/Documents/PhD/GenderAnalysis/TP53_interactions/Candidates.2018.03.01.csv")

p53.high = p53.int[which(p53.int$score>=0.6),'ID']
p53.med = p53.int[which(p53.int$score>=0.4),'ID']
p53.low = p53.int$ID

# p53.int.string = c("EFNB1","PSMD10","STAG2","MAGED2","EDA2R","RNF128","HUWE1","TAF1",
#             "BMX","MAGEB18","TSC22D3","MAGEA2B","ACSL4","MAGEH1","PLAC1","BRCC3",
#             "SH2D1A","UBE2A","SLC25A5","ATRX","RPS4X","AR","APEX2","FOXP3","RBM3",
#             "ARX","POLA1","PHKA2","UTP14A","DDX3X","CUL4B","CITED1","YY2")
# p53.int.boa = c("XGY2","PRKX","MIR4767","VCX","TBL1X","WWC3","NHS-AS1","NHS","PHKA2-AS1",
#                 "SH3KBP1","MAGEB5","TAB3","DMD","RNU6-16P","TMEM47","BCOR","NYX","PPP1R2P9",
#                 "LOC101927501","KDM6A","LOC401585","LINC01186","JADE3","ZNF81","PHF8","UQCRBP1",
#                 "SPIN4","MIR223","EDA2R","IGBP1","SNX12","CXCR3","RGAG4","ITM2A","TGIF2LX",
#                 "FAM133A","RPA4","DIAPH2-AS1","TSPAN6","ARL13A","NGFRAP1","TEX13B","ATG4A","COL4A6",
#                 "ACSL4","AMMECR1","HTR2C","PLS3","SLC25A5","XIAP","SMARCA1","GPC3","PLAC1","FAM122B",
#                 "FGF13","MIR320D2","SPANXB1","LINC00894","MIR4330","PASD1","GABRE","FLNA")

###Read DNA VAF#####
DNA.VAF = read.csv("~/Documents/PhD/GenderAnalysis/TCGA/Analysis/full.reduced.all.raw.TCGA.curated.mutations.csv", as.is=T)
DNA.VAF %>% filter(CANCER_TYPE%in%datasets) -> DNA.VAF
DNA.VAF %>% group_by(GENDER,CANCER_TYPE,PATIENT_ID) %>% dplyr::summarise() %>% 
  group_by(GENDER,CANCER_TYPE) %>% dplyr::summarise(N.P = n()) -> patient.numbers.by.cancer

DNA.VAF %>% group_by(GENDER,PATIENT_ID) %>% dplyr::summarise() %>% 
  group_by(GENDER) %>% dplyr::summarise(N.P = n()) -> patient.numbers


p53.mutants = filter(DNA.VAF,HUGO_SYMBOL=='TP53')
p53.mutants = unique(p53.mutants$PATIENT_ID)

filter.out = c("Silent",'Intron','IGR',"In_Frame_Ins" ,"In_Frame_Del", "lincRNA" )

filter(DNA.VAF,!(VARIANT_CLASSIFICATION%in%filter.out)) -> DNA.VAF

DNA.VAF$TP53.STATUS = 'Wt'
DNA.VAF[which(DNA.VAF$PATIENT_ID%in%p53.mutants),'TP53.STATUS'] = 'Mt'

DNA.VAF %>% group_by(GENDER,CANCER_TYPE,PATIENT_ID,TP53.STATUS) %>% dplyr::summarise() %>% 
  group_by(GENDER,CANCER_TYPE,TP53.STATUS) %>% dplyr::summarise(N.P = n()) -> patient.numbers.by.cancer.tp53.status

DNA.VAF %>% group_by(GENDER,PATIENT_ID,TP53.STATUS) %>% dplyr::summarise() %>% 
  group_by(GENDER,TP53.STATUS) %>% dplyr::summarise(N.P = n()) -> patient.numbers.tp53.status

#Use mutations present in 10+ patients and discard patients with too many mutations 10000+
DNA.VAF %>% group_by( PATIENT_ID) %>% dplyr::summarise(n=n()) %>% filter(n<=5000) -> aux
DNA.VAF %>% filter(PATIENT_ID %in% aux$PATIENT_ID) -> DNA.VAF

DNA.VAF %>% filter(CHROMOSOME=='X') -> DNA.X

#Use genes with 10+ mutations
DNA.X %>% group_by(HUGO_SYMBOL, PATIENT_ID) %>% dplyr::summarise(s=max(START_POSITION)) %>% 
  dplyr::mutate(U_ID = paste(HUGO_SYMBOL, PATIENT_ID,s,sep = '.')) -> aux

DNA.X$U_ID = paste(RMAF$HUGO_SYMBOL, RMAF$PATIENT_ID,RMAF$START_POSITION,sep = '.')
DNA.X %>% filter(U_ID %in% aux$U_ID) -> DNA.X

DNA.X %>% group_by(HUGO_SYMBOL) %>% dplyr::summarise(n=n()) %>% filter(n>=10) -> aux
DNA.X %>% filter(HUGO_SYMBOL %in% aux$HUGO_SYMBOL) -> RMAF
rm(aux)




################################################################################################
###############Build SM and EM Matrix for PS tests#############################
################################################################################################
badchars <- "[\xb5]|[\n]|[,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[/]|[\\]|[.]|[_]|[ ]"
clinical.data$PATHOLOGIC_STAGE <- toupper(gsub(pattern=badchars, replacement=".", x=clinical.data$PATHOLOGIC_STAGE))
PS.SM.mat = data.frame(matrix(data = 0,nrow = length(unique(DNA.VAF$PATIENT_ID)), 
                       ncol = length(unique(DNA.X$HUGO_SYMBOL))),
                       stringsAsFactors = F, row.names = unique(DNA.VAF$PATIENT_ID))
colnames(PS.SM.mat) = unique(DNA.X$HUGO_SYMBOL)
for (g in unique(DNA.X$HUGO_SYMBOL)){
  PS.SM.mat[which(rownames(PS.SM.mat)%in%(DNA.X[which(DNA.X$HUGO_SYMBOL==g),'PATIENT_ID'])),g]  = 1
  PS.SM.mat[which(rownames(PS.SM.mat)%in%(DNA.X[which(DNA.X$HUGO_SYMBOL==g&DNA.X$GENDER=='male'),'PATIENT_ID'])),g]  = 2
}

# PS.SM.mat = data.frame(matrix(data = NA,nrow = length(unique(DNA.VAF$PATIENT_ID)), 
#                               ncol = length(unique(RMAF$HUGO_SYMBOL))),
#                        stringsAsFactors = F, row.names = unique(DNA.VAF$PATIENT_ID))
# colnames(PS.SM.mat) = unique(RMAF$HUGO_SYMBOL)
# for (g in unique(RMAF$HUGO_SYMBOL)){
#   PS.SM.mat[which(rownames(PS.SM.mat)%in%(RMAF[which(RMAF$SM.RMAF==1&RMAF$HUGO_SYMBOL==g),'PATIENT_ID'])),g]  = 1
#   PS.SM.mat[which(rownames(PS.SM.mat)%in%(RMAF[which(RMAF$SM.RMAF==0&RMAF$HUGO_SYMBOL==g),'PATIENT_ID'])),g]  = 0
# }


#########################################################
#####Do Statistical test with Propensity scores##########
#########################################################
cf = cofounding.factors.ind
cfd = cofounding.factors.default
####DO for SM
annot = clinical.data[which(clinical.data$PATIENT_ID%in%rownames(PS.SM.mat)),unique(c('PATIENT_ID',cf,cfd))]
keep = rownames(PS.SM.mat)%in%clinical.data$PATIENT_ID
PS.SM.mat = PS.SM.mat[keep,]
res = DoStatTest(clinical.data = annot,molecular.data = PS.SM.mat,default.confounding.factors = cfd, 
                 extra.cofounding.factors = c('SMOKING_STATUS','RACE'), is.continuous = T)
res$DATASET = 'ALL'
res = as.data.frame(res)
#Add SM numbers and total numbers 
tot.fem = length(rownames(PS.SM.mat)[which(rownames(PS.SM.mat)%in%filter(clinical.data,GENDER=='female')$PATIENT_ID)])
tot.male = length(rownames(PS.SM.mat)[which(rownames(PS.SM.mat)%in%filter(clinical.data,GENDER=='male')$PATIENT_ID)])
res$total.samples = paste(tot.male,'M',":","F",tot.fem)
res$total.males = tot.male
res$total.females = tot.fem
res$mutations = ''
res$mutations.males = 0
res$mutations.females = 0
for(g in colnames(PS.SM.mat)){
  sm.fem = sum(PS.SM.mat[which(rownames(PS.SM.mat)%in%filter(clinical.data,GENDER=='female')$PATIENT_ID),g],na.rm = T)
  sm.male = sum(PS.SM.mat[which(rownames(PS.SM.mat)%in%filter(clinical.data,GENDER=='male')$PATIENT_ID),g],na.rm = T)
  res[which(res$feature==g),'mutations'] = paste(sm.male,'M',":","F",sm.fem)
  res[which(res$feature==g),'mutations.males'] = sm.male
  res[which(res$feature==g),'mutations.females'] = sm.fem
  
}

for (cancer in datasets){
  print(cancer)
  DNA.X %>% filter(CANCER_TYPE==cancer) -> RMAF.SM
  DNA.VAF %>% filter(CANCER_TYPE==cancer) -> DNA.cancer
  PS.SM.mat = data.frame(matrix(data = 0,nrow = length(unique(DNA.cancer$PATIENT_ID)), 
                                ncol = length(unique(RMAF.SM$HUGO_SYMBOL))),
                         stringsAsFactors = F, row.names = unique(DNA.cancer$PATIENT_ID))
  colnames(PS.SM.mat) = unique(RMAF.SM$HUGO_SYMBOL)
  for (g in unique(RMAF.SM$HUGO_SYMBOL)){
    RMAF.SM %>% filter(HUGO_SYMBOL==g) -> RMAF.g
    PS.SM.mat[which(rownames(PS.SM.mat)%in%RMAF.g$PATIENT_ID),g] = 1
    PS.SM.mat[which(rownames(PS.SM.mat)%in%RMAF.g[which(RMAF.g$GENDER=='male'),'PATIENT_ID']),g] = 2
  }
  clinical.data %>% filter(CANCER_TYPE==cancer) -> clinical.cancer
  annot = clinical.cancer[which(clinical.cancer$PATIENT_ID%in%rownames(PS.SM.mat)),unique(c('PATIENT_ID',cf,cfd))]
  keep = rownames(PS.SM.mat)%in%clinical.cancer$PATIENT_ID
  PS.SM.mat = PS.SM.mat[keep,]
  res.aux = DoStatTest(clinical.data = annot,molecular.data = PS.SM.mat, is.continuous = T,
                       default.confounding.factors = cfd, extra.cofounding.factors = cf)
  res.aux = as.data.frame(res.aux)
  res.aux$DATASET = cancer
  ####ADD TOTAL NUMBERS
  tot.fem = length(rownames(PS.SM.mat)[which(rownames(PS.SM.mat)%in%filter(clinical.data,GENDER=='female')$PATIENT_ID)])
  tot.male = length(rownames(PS.SM.mat)[which(rownames(PS.SM.mat)%in%filter(clinical.data,GENDER=='male')$PATIENT_ID)])
  res.aux$total.samples = paste(tot.male,'M',":","F",tot.fem)
  res.aux$total.males = tot.male
  res.aux$total.females = tot.fem
  res.aux$mutations.males = 0
  res.aux$mutations.females = 0
  res.aux$mutations = ''
  for(g in colnames(PS.SM.mat)){
    sm.fem = sum(PS.SM.mat[which(rownames(PS.SM.mat)%in%filter(clinical.data,GENDER=='female')$PATIENT_ID),g])
    sm.male = sum(PS.SM.mat[which(rownames(PS.SM.mat)%in%filter(clinical.data,GENDER=='male')$PATIENT_ID),g])
    res.aux[which(res.aux$feature==g),'mutations'] = paste(sm.male,'M',":","F",sm.fem)
    res.aux[which(res.aux$feature==g),'mutations.males'] = sm.male
    res.aux[which(res.aux$feature==g),'mutations.females'] = sm.fem
  }
  
  ###Merge data.frame
  res = rbind(res,as.data.frame(res.aux))
}

write.csv(res,paste(output.dir,'AllDataSets/','DNA.adj.muts.Chi.Test.Propensity.csv',sep = ''),row.names = F)


###Do p53 interactors tests
res$z.score = sign(res$coef)*sqrt(res$chisq)
res$TP53.INT.HIGH = 'No'
res[which(res$feature%in%p53.high),'TP53.INT.HIGH'] = 'Yes'
res$TP53.INT.MED = 'No'
res[which(res$feature%in%p53.med),'TP53.INT.MED'] = 'Yes'
res$TP53.INT.LOW = 'No'
res[which(res$feature%in%p53.low),'TP53.INT.LOW'] = 'Yes'
test.tbl = NULL


for (cancer in c('ALL',datasets)){
  tmp = filter(res,DATASET==cancer)
  test.high = t.test(tmp$z.score~tmp$TP53.INT.HIGH)
  test.med =  t.test(tmp$z.score~tmp$TP53.INT.MED)
  test.low =  t.test(tmp$z.score~tmp$TP53.INT.LOW)
  if (is.null(test.tbl)){
    test.tbl = data.frame(cancer,test.high$statistic,test.high$p.value,'High.conf')
    colnames(test.tbl) = c("DATASET",'t-statistic','p-value','Interactors')
  }
  else{
    test.tmp = data.frame(cancer,test.high$statistic,test.high$p.value,'High.conf')
    colnames(test.tmp) = c("DATASET",'t-statistic','p-value','Interactors')
    test.tbl = rbind(test.tbl,test.tmp)
  }
  
  test.tmp = data.frame(cancer,test.med$statistic,test.med$p.value,'Med.conf')
  colnames(test.tmp) = c("DATASET",'t-statistic','p-value','Interactors')
  test.tbl = rbind(test.tbl,test.tmp)
  test.tmp = data.frame(cancer,test.low$statistic,test.low$p.value,'Low.conf')
  colnames(test.tmp) = c("DATASET",'t-statistic','p-value','Interactors')
  test.tbl = rbind(test.tbl,test.tmp)
  
}
write.csv(test.tbl,paste(output.dir,'AllDataSets/','p53.int.DNA.muts.t.test.csv',sep = ''),row.names = F)


####Plot Density of Chi-squares ##############
DoDensityPlotFacetGenderRMAF(dat = res, dat.col = 'z.score', facet.col = 'DATASET'
                             ,fill.col = 'TP53.INT.HIGH',
                             title.txt = paste("TP53 Status SM Z-Scores Densities"))  
DoDensityPlotFacetGenderRMAF(dat = res, dat.col = 'z.score', facet.col = 'DATASET'
                             ,fill.col = 'TP53.INT.MED',
                             title.txt = paste("TP53 Status SM Z-Scores Densities"))  
DoDensityPlotFacetGenderRMAF(dat = res, dat.col = 'z.score', facet.col = 'DATASET'
                             ,fill.col = 'TP53.INT.LOW',
                             title.txt = paste("TP53 Status SM Z-Scores Densities"))  


DoDensityPlotGenderRMAF(dat = filter(res,DATASET=='ALL'),dat.col = 'z.score',fill.col = 'TP53.INT.HIGH',title.txt = 'epa')
DoDensityPlotGenderRMAF(dat = filter(res,DATASET=='ALL'),dat.col = 'z.score',fill.col = 'TP53.INT.LOW',title.txt = 'epa')

##############################################

###Do for EM
RMAF %>% filter(EM.RMAF==1) -> RMAF.EM
PS.EM.mat = data.frame(matrix(data = 0,nrow = length(unique(DNA.VAF$PATIENT_ID)), 
                              ncol = length(unique(RMAF.EM$HUGO_SYMBOL))),
                       stringsAsFactors = F, row.names = unique(DNA.VAF$PATIENT_ID))
colnames(PS.EM.mat) = unique(RMAF.EM$HUGO_SYMBOL)
for (g in unique(RMAF.EM$HUGO_SYMBOL)){
  RMAF.EM %>% filter(HUGO_SYMBOL==g) -> RMAF.g
  PS.EM.mat[which(rownames(PS.EM.mat)%in%RMAF.g$PATIENT_ID),g] = 1
}

annot = clinical.data[which(clinical.data$PATIENT_ID%in%rownames(PS.EM.mat)),unique(c('PATIENT_ID',cf,cfd))]
keep = rownames(PS.EM.mat)%in%clinical.data$PATIENT_ID
PS.EM.mat = PS.EM.mat[keep,]
res = DoStatTest(clinical.data = annot,molecular.data = PS.EM.mat,default.confounding.factors = cfd, extra.cofounding.factors = c('SMOKING_STATUS','RACE'))
res$DATASET = 'ALL'
res = as.data.frame(res)
#Add EM numbers and total numbers 
tot.fem = length(rownames(PS.EM.mat)[which(rownames(PS.EM.mat)%in%filter(clinical.data,GENDER=='female')$PATIENT_ID)])
tot.male = length(rownames(PS.EM.mat)[which(rownames(PS.EM.mat)%in%filter(clinical.data,GENDER=='male')$PATIENT_ID)])
res$total.samples = paste(tot.male,'M',":","F",tot.fem)
res$total.females = tot.fem
res$total.males = tot.male
res$expressed.mutations = ''
res$expressed.mutations.females = 0
res$expressed.mutations.males = 0

for(g in colnames(PS.EM.mat)){
  em.fem = sum(PS.EM.mat[which(rownames(PS.EM.mat)%in%filter(clinical.data,GENDER=='female')$PATIENT_ID),g])
  em.male = sum(PS.EM.mat[which(rownames(PS.EM.mat)%in%filter(clinical.data,GENDER=='male')$PATIENT_ID),g])
  res[which(res$feature==g),'expressed.mutations'] = paste(em.male,'M',":","F",em.fem)
  res[which(res$feature==g),'expressed.mutations.females'] = em.fem
  res[which(res$feature==g),'expressed.mutations.males'] = em.male
  
}
for (cancer in datasets){
  print(cancer)
  RMAF %>% filter(EM.RMAF==1,CANCER_TYPE==cancer) -> RMAF.EM
  DNA.VAF %>% filter(CANCER_TYPE==cancer) -> DNA.cancer
  PS.EM.mat = data.frame(matrix(data = 0,nrow = length(unique(DNA.cancer$PATIENT_ID)), 
                                ncol = length(unique(RMAF.EM$HUGO_SYMBOL))),
                         stringsAsFactors = F, row.names = unique(DNA.cancer$PATIENT_ID))
  colnames(PS.EM.mat) = unique(RMAF.EM$HUGO_SYMBOL)
  for (g in unique(RMAF.EM$HUGO_SYMBOL)){
    RMAF.EM %>% filter(HUGO_SYMBOL==g) -> RMAF.g
    PS.EM.mat[which(rownames(PS.EM.mat)%in%RMAF.g$PATIENT_ID),g] = 1
  }
  clinical.data %>% filter(CANCER_TYPE==cancer) -> clinical.cancer
  annot = clinical.cancer[which(clinical.cancer$PATIENT_ID%in%rownames(PS.EM.mat)),unique(c('PATIENT_ID',cf,cfd))]
  keep = rownames(PS.EM.mat)%in%clinical.cancer$PATIENT_ID
  PS.EM.mat = PS.EM.mat[keep,]
  res.aux = DoStatTest(clinical.data = annot,molecular.data = PS.EM.mat,default.confounding.factors = cfd, extra.cofounding.factors = cf)
  res.aux$DATASET = cancer
  res.aux = as.data.frame(res.aux)
  ####ADD TOTAL NUMBERS
  tot.fem = length(rownames(PS.EM.mat)[which(rownames(PS.EM.mat)%in%filter(clinical.data,GENDER=='female')$PATIENT_ID)])
  tot.male = length(rownames(PS.EM.mat)[which(rownames(PS.EM.mat)%in%filter(clinical.data,GENDER=='male')$PATIENT_ID)])
  res.aux$total.samples = paste(tot.male,'M',":","F",tot.fem)
  res.aux$total.females = tot.fem
  res.aux$total.males = tot.male
  res.aux$expressed.mutations = ''
  res.aux$expressed.mutations.females = 0
  res.aux$expressed.mutations.males = 0
  for(g in colnames(PS.EM.mat)){
    sm.fem = sum(PS.EM.mat[which(rownames(PS.EM.mat)%in%filter(clinical.data,GENDER=='female')$PATIENT_ID),g])
    sm.male = sum(PS.EM.mat[which(rownames(PS.EM.mat)%in%filter(clinical.data,GENDER=='male')$PATIENT_ID),g])
    res.aux[which(res.aux$feature==g),'expressed.mutations'] = paste(sm.male,'M',":","F",sm.fem)
    res.aux[which(res.aux$feature==g),'expressed.mutations.females'] = em.fem
    res.aux[which(res.aux$feature==g),'expressed.mutations.males'] = em.male
  }
  
  ###Merge data.frame
  res = rbind(res,as.data.frame(res.aux))
  
}

write.csv(res,paste(output.dir,'AllDataSets/','EM.Chi.Test.Propensity.csv',sep = ''),row.names = F)


###Do p53 interactors tests
res$z.score = sign(res$coef)*sqrt(res$chisq)
res$TP53.INT.HIGH = 'No'
res[which(res$feature%in%p53.high),'TP53.INT.HIGH'] = 'Yes'
res$TP53.INT.MED = 'No'
res[which(res$feature%in%p53.med),'TP53.INT.MED'] = 'Yes'
res$TP53.INT.LOW = 'No'
res[which(res$feature%in%p53.low),'TP53.INT.LOW'] = 'Yes'
test.tbl = NULL


for (cancer in c('ALL')){
  tmp = filter(res,DATASET==cancer)
  test.high = t.test(tmp$z.score~tmp$TP53.INT.HIGH)
  test.med =  t.test(tmp$z.score~tmp$TP53.INT.MED)
  test.low =  t.test(tmp$z.score~tmp$TP53.INT.LOW)
  if (is.null(test.tbl)){
    test.tbl = data.frame(cancer,test.high$statistic,test.high$p.value,'High.conf')
    colnames(test.tbl) = c("DATASET",'t-statistic','p-value','Interactors')
  }
  else{
    test.tmp = data.frame(cancer,test.high$statistic,test.high$p.value,'High.conf')
    colnames(test.tmp) = c("DATASET",'t-statistic','p-value','Interactors')
    test.tbl = rbind(test.tbl,test.tmp)
  }
  
  test.tmp = data.frame(cancer,test.med$statistic,test.med$p.value,'Med.conf')
  colnames(test.tmp) = c("DATASET",'t-statistic','p-value','Interactors')
  test.tbl = rbind(test.tbl,test.tmp)
  test.tmp = data.frame(cancer,test.low$statistic,test.low$p.value,'Low.conf')
  colnames(test.tmp) = c("DATASET",'t-statistic','p-value','Interactors')
  test.tbl = rbind(test.tbl,test.tmp)
  
}
write.csv(test.tbl,paste(output.dir,'AllDataSets/','p53.int.EM.t.test.csv',sep = ''),row.names = F)



#############TEST MT-p53 EFFECT###################
RMAF %>% filter(SM.RMAF==1,GENDER=='female') -> RMAF.SM
PS.SM.mat = data.frame(matrix(data = 0,nrow = length(unique(DNA.VAF[which(DNA.VAF$GENDER=='female'),"PATIENT_ID"])), 
                              ncol = length(unique(RMAF.SM$HUGO_SYMBOL))),
                       stringsAsFactors = F, row.names = unique(DNA.VAF[which(DNA.VAF$GENDER=='female'),"PATIENT_ID"]))
colnames(PS.SM.mat) = unique(RMAF.SM$HUGO_SYMBOL)
for (g in unique(RMAF.SM$HUGO_SYMBOL)){
  RMAF.SM %>% filter(HUGO_SYMBOL==g) -> RMAF.g
  PS.SM.mat[which(rownames(PS.SM.mat)%in%RMAF.g$PATIENT_ID),g] = 1
}

annot = clinical.data[which(clinical.data$PATIENT_ID%in%rownames(PS.SM.mat)),unique(c('PATIENT_ID',cf,cfd))]
keep = rownames(PS.SM.mat)%in%clinical.data$PATIENT_ID
PS.SM.mat = PS.SM.mat[keep,]
###Hardcode for Mt-p53
annot$GENDER = ifelse(annot$PATIENT_ID%in%p53.mutants,'female','male')
res = DoStatTest(clinical.data = annot,molecular.data = PS.SM.mat,default.confounding.factors = cfd, extra.cofounding.factors = c('SMOKING_STATUS','RACE'))
res = as.data.frame(res)
res$z.score = sign(res$coef)*sqrt(res$chisq)
#Add SM numbers and total numbers 
clinical.data$TP53_STATUS = 'Wt'
clinical.data[which(clinical.data$PATIENT_ID%in%p53.mutants),"TP53_STATUS"] = 'Mt'

tot.mt.fem = length(rownames(PS.SM.mat)[which(rownames(PS.SM.mat)%in%filter(clinical.data,GENDER=='female',TP53_STATUS=='Mt')$PATIENT_ID)])
tot.wt.fem = length(rownames(PS.SM.mat)[which(rownames(PS.SM.mat)%in%filter(clinical.data,GENDER=='female',TP53_STATUS=='Wt')$PATIENT_ID)])
res$total.samples = paste(tot.wt.fem,'F.Wt-p53',":","F.Mt-p53",tot.mt.fem)
res$silent.mutations = ''
for(g in colnames(PS.SM.mat)){
  sm.mt.fem = sum(PS.SM.mat[which(rownames(PS.SM.mat)%in%filter(clinical.data,GENDER=='female',TP53_STATUS=='Mt')$PATIENT_ID),g])
  sm.wt.fem = sum(PS.SM.mat[which(rownames(PS.SM.mat)%in%filter(clinical.data,GENDER=='female',TP53_STATUS=='Wt')$PATIENT_ID),g])
  res[which(res$feature==g),'silent.mutations'] = paste(sm.wt.fem,'F.Wt-p53',":","F.Mt-p53",sm.mt.fem)
}
write.csv(res,paste(output.dir,'AllDataSets/','p53.status.females.SM.chi.sq.test.csv',sep = ''),row.names = F)

plot(density(rnorm(n = 531,mean = 0,sd = 1)),main = 'Black: Normal Dist (mean = 0, sd = 1, n = 531) \n Red: Mt-p53 vs Wt-p53 FEMALES test \n(mean = 0.47 , sd = 0.99, n = 531) ')
lines(density(res$z.score),col='red')
lines(density(res[which(res$feature%in%p53.int),'z.score']),col='blue')

plot(density(res[which(!(res$feature%in%p53.int)),'z.score']),col='blue')
lines(density(res[which(res$feature%in%p53.int),'z.score']),col='red')


##################################################


################################################################################################
###############Silenced Mutation Tests#############################
################################################################################################
#patient.numbers$GENDER = factor(patient.numbers$GENDER,levels=c("female","male"))
# res.SM = NULL
# for (g in unique(RMAF$HUGO_SYMBOL)){
#   
#   #print(g)
#   #DO exact test
#   RMAF %>% filter(HUGO_SYMBOL==g,SM.RMAF==1) -> RMAF.g
#   RMAF.g$GENDER = factor(RMAF.g$GENDER,levels=c("female","male"))
#   
#   mat = data.frame(N.EM=table(RMAF.g$GENDER[!duplicated(RMAF.g$PATIENT_ID)]), N.P=patient.numbers$N.P)
#   colnames(mat) = c('GENDER','NP.EM','NP')
#   
#   RMAF.g %>% group_by(GENDER) %>% dplyr::summarise(AVG.RMAF=mean(RMAF)) -> RMAF.avg
#      
#   # 
#   #patient.numbers %>% left_join(numbers.g) -> mat
#   
#   mat$NP.NO.EM = mat$NP - mat$NP.EM
#   t = prop.test(mat$NP.EM,mat$NP)
#   prop.stat = t$statistic
#   prop.pvalue = t$p.value
#   mat$NP = NULL
#   #print(mat)
#  
#   p.val = fisher.test(mat[,2:3],alternative = 'two.sided')
#   z.stat = Exact_Statistic(mat[1,2],mat[2,2],mat[1,3],mat[2,3])
#   
#   mat %>% left_join(RMAF.avg) -> mat.aux
#   res.aux = data.frame(g,p.val$p.value,z.stat,as.double(prop.stat),as.double(prop.pvalue),mat[1,2]/mat[1,3]*100,mat[2,2]/mat[2,3]*100,mat[1,2],mat[2,2],
#                        mat.aux[which(mat.aux$GENDER=='female'),'AVG.RMAF'],
#                        mat.aux[which(mat.aux$GENDER=='male'),'AVG.RMAF'],
#                        'ALL')
#   colnames(res.aux) = c("HUGO_SYMBOL","P.VALUE.FISHER",'Z.STAT',"PROP.STAT","PROP.P.VALUE","FEMALE.SM.PRCT","MALE.SM.PRCT",
#                        "FEMALES.W.SM","MALES.W.SM","FEMALE.AVG.RMAF","MALE.AVG.RMAF","DATASET")
#   
#   if(is.null(res.SM)){
#     res.SM <- res.aux
#   }
#   else{
#     res.SM <- rbind(res.SM,res.aux)  
#   }
#   
# }
# #write.csv(res,paste(output.dir,'AllDataSets/','ExactTestsSM.csv',sep = ''),row.names = F)
# res.SM$TP53.INT = 0
# res.SM[which(res.SM$HUGO_SYMBOL%in%p53.int),'TP53.INT'] = 1
# test = t.test(res.SM$Z.STAT~res.SM$TP53.INT)
# print('P.VAL All cancers SM')
# print(test$p.value)
# 
# ##############################################################################################
# ################################################################################################
# ###############Expressed Mutation Tests#############################
# ################################################################################################
# ################################################################################################
# #patient.numbers$GENDER = factor(patient.numbers$GENDER,levels=c("female","male"))
# res.EM= NULL
# for (g in unique(RMAF$HUGO_SYMBOL)){
#   
#   #print(g)
#   #DO exact test
#   RMAF %>% filter(HUGO_SYMBOL==g,EM.RMAF==1) -> RMAF.g
#   RMAF.g$GENDER = factor(RMAF.g$GENDER,levels=c("female","male"))
#   
#   mat = data.frame(N.EM=table(RMAF.g$GENDER[!duplicated(RMAF.g$PATIENT_ID)]), N.P=patient.numbers$N.P)
#   colnames(mat) = c('GENDER','NP.EM','NP')
#   
#   RMAF.g %>% group_by(GENDER) %>% dplyr::summarise(AVG.RMAF=mean(RMAF)) -> RMAF.avg
#   
#   # 
#   #patient.numbers %>% left_join(numbers.g) -> mat
#   
#   mat$NP.NO.EM = mat$NP - mat$NP.EM
#   mat$NP = NULL
#   #print(mat)
#   
#   p.val.fish = fisher.test(mat[,2:3],alternative = 'two.sided')
#  
#   mat %>% left_join(RMAF.avg) -> mat.aux
#   res.aux = data.frame(g,p.val$p.value,mat[1,2]/mat[1,3]*100,mat[2,2]/mat[2,3]*100,mat[1,2],mat[2,2],
#                        mat.aux[which(mat.aux$GENDER=='female'),'AVG.RMAF'],
#                        mat.aux[which(mat.aux$GENDER=='male'),'AVG.RMAF'],
#                        'ALL')
#   colnames(res.aux) = c("HUGO_SYMBOL","P.VALUE","FEMALE.EM.PRCT","MALE.EM.PRCT",
#                         "FEMALES.W.EM","MALES.W.EM","FEMALE.AVG.RMAF","MALE.AVG.RMAF","DATASET")
#   
#   if(is.null(res.EM)){
#     res.EM <- res.aux
#   }
#   else{
#     res.EM <- rbind(res.EM,res.aux)  
#   }
#   
# }
# 
# ################################################################################################
# ###############Silenced Mutation Tests#############################
# ################################################################################################
# #patient.numbers$GENDER = factor(patient.numbers$GENDER,levels=c("female","male"))
# for (cancer in datasets){
#   for (g in unique(RMAF$HUGO_SYMBOL)){
#     
#     #print(g)
#     #DO exact test
#     RMAF %>% filter(HUGO_SYMBOL==g,SM.RMAF==1, CANCER_TYPE==cancer) -> RMAF.g
#     RMAF.g$GENDER = factor(RMAF.g$GENDER,levels=c("female","male"))
#     
#     mat = data.frame(N.EM=table(RMAF.g$GENDER[!duplicated(RMAF.g$PATIENT_ID)]), N.P=patient.numbers$N.P)
#     colnames(mat) = c('GENDER','NP.EM','NP')
#     
#     RMAF.g %>% group_by(GENDER) %>% dplyr::summarise(AVG.RMAF=mean(RMAF)) -> RMAF.avg
#     
#     # 
#     #patient.numbers %>% left_join(numbers.g) -> mat
#     
#     mat$NP.NO.EM = mat$NP - mat$NP.EM
#     mat$NP = NULL
#     #print(mat)
#     
#     p.val = fisher.test(mat[,2:3],alternative = 'two.sided')
#     mat %>% left_join(RMAF.avg) -> mat.aux
#     res.aux = data.frame(g,p.val$p.value,mat[1,2]/mat[1,3]*100,mat[2,2]/mat[2,3]*100,mat[1,2],mat[2,2],
#                          mat.aux[which(mat.aux$GENDER=='female'),'AVG.RMAF'],
#                          mat.aux[which(mat.aux$GENDER=='male'),'AVG.RMAF'],
#                          cancer)
#     colnames(res.aux) = c("HUGO_SYMBOL","P.VALUE","FEMALE.SM.PRCT","MALE.SM.PRCT",
#                           "FEMALES.W.SM","MALES.W.SM","FEMALE.AVG.RMAF","MALE.AVG.RMAF","DATASET")
#     
#     if(is.null(res.SM)){
#       res.SM <- res.aux
#     }
#     else{
#       res.SM <- rbind(res.SM,res.aux)  
#     }
#     
#   }
# }
# 
# #write.csv(res,paste(output.dir,'AllDataSets/','ExactTestsSM.csv',sep = ''),row.names = F)
# # res.SM$TP53.INT = 0
# # res.SM[which(res.SM$HUGO_SYMBOL%in%p53.int),'TP53.INT'] = 1
# # test = t.test(res.SM$P.VALUE~res.SM$TP53.INT)
# # print('P.VAL All cancers SM')
# # print(test$p.value)
# 
# ##############################################################################################
# ################################################################################################
# ###############Expressed Mutation Tests#############################
# ################################################################################################
# ################################################################################################
# #patient.numbers$GENDER = factor(patient.numbers$GENDER,levels=c("female","male"))
# for (cancer in datasets){
#   for (g in unique(RMAF$HUGO_SYMBOL)){
#     
#     #print(g)
#     #DO exact test
#     RMAF %>% filter(HUGO_SYMBOL==g,EM.RMAF==1, CANCER_TYPE==cancer) -> RMAF.g
#     RMAF.g$GENDER = factor(RMAF.g$GENDER,levels=c("female","male"))
#     
#     mat = data.frame(N.EM=table(RMAF.g$GENDER[!duplicated(RMAF.g$PATIENT_ID)]), N.P=patient.numbers$N.P)
#     colnames(mat) = c('GENDER','NP.EM','NP')
#     
#     RMAF.g %>% group_by(GENDER) %>% dplyr::summarise(AVG.RMAF=mean(RMAF)) -> RMAF.avg
#     
#     # 
#     #patient.numbers %>% left_join(numbers.g) -> mat
#     
#     mat$NP.NO.EM = mat$NP - mat$NP.EM
#     mat$NP = NULL
#     #print(mat)
#     
#     p.val = fisher.test(mat[,2:3],alternative = 'two.sided')
#     mat %>% left_join(RMAF.avg) -> mat.aux
#     res.aux = data.frame(g,p.val$p.value,mat[1,2]/mat[1,3]*100,mat[2,2]/mat[2,3]*100,mat[1,2],mat[2,2],
#                          mat.aux[which(mat.aux$GENDER=='female'),'AVG.RMAF'],
#                          mat.aux[which(mat.aux$GENDER=='male'),'AVG.RMAF'],
#                          cancer)
#     colnames(res.aux) = c("HUGO_SYMBOL","P.VALUE","FEMALE.EM.PRCT","MALE.EM.PRCT",
#                           "FEMALES.W.EM","MALES.W.EM","FEMALE.AVG.RMAF","MALE.AVG.RMAF","DATASET")
#     
#     if(is.null(res.EM)){
#       res.EM <- res.aux
#     }
#     else{
#       res.EM <- rbind(res.EM,res.aux)  
#     }
#     
#   }
# }
# 
# write.csv(res.SM,paste(output.dir,'AllDataSets/','ExactTestsSM.csv',sep = ''),row.names = F)
# write.csv(res.EM,paste(output.dir,'AllDataSets/','ExactTestsEM.csv',sep = ''),row.names = F)
# 
# write.csv(res.SM[which(res.SM$HUGO_SYMBOL%in%p53.int),],paste(output.dir,'AllDataSets/','ExactTestsSMp53Int.csv',sep = ''),row.names = F)
# write.csv(res.EM[which(res.EM$HUGO_SYMBOL%in%p53.int),],paste(output.dir,'AllDataSets/','ExactTestsEMp53Int.csv',sep = ''),row.names = F)

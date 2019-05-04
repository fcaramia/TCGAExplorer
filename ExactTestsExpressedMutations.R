rm(list = ls())
library(dplyr)
library(gtools)
library(data.table)
source('PropensityScoresFunctions.R')
source('ExpressionPlotsFunctions.R')
Exact_Statistic <- function(a, b, c , d){
  
  n = (a+1/2)*(c+1/2)/(b+1/2)*(d+1/2)
  d = sqrt(1/(a+1/2)+1/(b+1/2)+1/(c+1/2)+1/(d+1/2))
  z = log(n)/d
  
  return(z)
}
#Datasets and directories and variables
datasets = c("ESCA","HNSC","LUSC","BLCA","LIHC","STAD","LGG","COAD","PAAD","READ","SKCM",'LUAD')
#datasets = c("KIRC")
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

p53.high = as.vector(p53.int[which(p53.int$score>=0.6),'ID'])
p53.med = as.vector(p53.int[which(p53.int$score>=0.4),'ID'])
p53.low = as.vector(p53.int$ID)

# p53.int.string = c("EFNB1","PSMD10","STAG2","MAGED2","EDA2R","RNF128","HUWE1","TAF1",
#              "BMX","MAGEB18","TSC22D3","MAGEA2B","ACSL4","MAGEH1","PLAC1","BRCC3",
#              "SH2D1A","UBE2A","SLC25A5","ATRX","RPS4X","AR","APEX2","FOXP3","RBM3",
#              "ARX","POLA1","PHKA2","UTP14A","DDX3X","CUL4B","CITED1","YY2")
# p53.int.boa = c("XGY2","PRKX","MIR4767","VCX","TBL1X","WWC3","NHS-AS1","NHS","PHKA2-AS1",
#                 "SH3KBP1","MAGEB5","TAB3","DMD","RNU6-16P","TMEM47","BCOR","NYX","PPP1R2P9",
#                 "LOC101927501","KDM6A","LOC401585","LINC01186","JADE3","ZNF81","PHF8","UQCRBP1",
#                 "SPIN4","MIR223","EDA2R","IGBP1","SNX12","CXCR3","RGAG4","ITM2A","TGIF2LX",
#                 "FAM133A","RPA4","DIAPH2-AS1","TSPAN6","ARL13A","NGFRAP1","TEX13B","ATG4A","COL4A6",
#                 "ACSL4","AMMECR1","HTR2C","PLS3","SLC25A5","XIAP","SMARCA1","GPC3","PLAC1","FAM122B",
#                 "FGF13","MIR320D2","SPANXB1","LINC00894","MIR4330","PASD1","GABRE","FLNA")

###Read DNA VAF#####
DNA.VAF = fread("~/Documents/PhD/GenderAnalysis/TCGA/Analysis/full.reduced.all.raw.TCGA.curated.mutations.csv")
DNA.VAF %>% filter(CANCER_TYPE%in%datasets) -> DNA.VAF
DNA.VAF %>% group_by(GENDER,CANCER_TYPE,PATIENT_ID) %>% dplyr::summarise() %>% 
  group_by(GENDER,CANCER_TYPE) %>% dplyr::summarise(N.P = n()) -> patient.numbers.by.cancer

write.csv(patient.numbers.by.cancer,paste(output.dir,'patient.numbers.by.cancer.csv',sep = ''),row.names = F)

DNA.VAF %>% group_by(GENDER,PATIENT_ID) %>% dplyr::summarise() %>% 
  group_by(GENDER) %>% dplyr::summarise(N.P = n()) -> patient.numbers


p53.mutants = filter(DNA.VAF,HUGO_SYMBOL=='TP53')
p53.mutants = unique(p53.mutants$PATIENT_ID)

filter.out = c("Silent",'Intron','IGR',"In_Frame_Ins" ,"In_Frame_Del", "lincRNA" )

filter(DNA.VAF,!(VARIANT_CLASSIFICATION%in%filter.out)) -> DNA.VAF

cancer.drivers = c('TP53','MYC','KRAS','PIK3CA')

DNA.VAF %>% filter(HUGO_SYMBOL%in%cancer.drivers) %>% group_by(GENDER,CANCER_TYPE,HUGO_SYMBOL,PATIENT_ID) %>% 
  dplyr::summarise()  %>%
  group_by(GENDER,CANCER_TYPE,HUGO_SYMBOL) %>% dplyr::summarise(N.M = n()) -> cancer.drivers.mutations

write.csv(cancer.drivers.mutations,paste(output.dir,'cancer.drivers.by.cancer.csv',sep = ''),row.names = F)

DNA.VAF  %>% group_by(GENDER,CANCER_TYPE,HUGO_SYMBOL,PATIENT_ID,CHROMOSOME) %>% 
  dplyr::summarise()  %>%
  group_by(GENDER,CANCER_TYPE,HUGO_SYMBOL,CHROMOSOME) %>% dplyr::summarise(N.M = n()) -> cancer.mutations.by.gene

cancer.mutations.by.gene$EXCEL = paste('`',cancer.mutations.by.gene$HUGO_SYMBOL,sep = '')

write.csv(cancer.mutations.by.gene,paste(output.dir,'mutation.by.gene.KIRC.csv',sep = ''),row.names = F)

DNA.VAF  %>% group_by(GENDER,HUGO_SYMBOL,PATIENT_ID,CHROMOSOME) %>% 
  dplyr::summarise()  %>%
  group_by(GENDER,HUGO_SYMBOL,CHROMOSOME) %>% dplyr::summarise(N.M = n()) -> cancer.mutations.by.gene.total
cancer.mutations.by.gene.total$EXCEL = paste('`',cancer.mutations.by.gene.total$HUGO_SYMBOL,sep='')
write.csv(cancer.mutations.by.gene.total,paste(output.dir,'mutation.by.gene.total.csv',sep = ''),row.names = F)


DNA.VAF$TP53.STATUS = 'Wt'
DNA.VAF[which(DNA.VAF$PATIENT_ID%in%p53.mutants),'TP53.STATUS'] = 'Mt'

DNA.VAF %>% group_by(GENDER,CANCER_TYPE,PATIENT_ID,TP53.STATUS) %>% dplyr::summarise() %>% 
  group_by(GENDER,CANCER_TYPE,TP53.STATUS) %>% dplyr::summarise(N.P = n()) -> patient.numbers.by.cancer.tp53.status

DNA.VAF %>% group_by(GENDER,PATIENT_ID,TP53.STATUS) %>% dplyr::summarise() %>% 
  group_by(GENDER,TP53.STATUS) %>% dplyr::summarise(N.P = n()) -> patient.numbers.tp53.status


######RMAF DATA#############
RMAF = read.csv("~/Documents/PhD/GenderAnalysis/TCGA/Analysis/all.TCGA.ExMut.V5.csv", as.is=T)

filter.out = c("Silent",'Intron','IGR',"In_Frame_Ins" ,"In_Frame_Del", "lincRNA" )

filter(RMAF,!is.na(ALT_COUNT),!(VARIANT_CLASSIFICATION%in%filter.out)) -> RMAF
filter(RMAF,!is.na(ALT_COUNT)) -> RMAF
#RMAF%>% filter(VARIANT_TYPE=='SNP') -> RMAF
RMAF %>% filter(CANCER_TYPE%in%datasets) -> RMAF
RMAF$TP53_STATUS = ifelse(RMAF$PATIENT_ID%in%p53.mutants,'Mt','Wt')
RMAF$TP53_INT = ifelse(RMAF$HUGO_SYMBOL%in%p53.int,'I','N')

RMAF$EM.RNA = 0
RMAF$EXPRS.GENE = 'EXPRS'
#Filter to only use cancer tissues

#Define Expressed mutation RNA
RMAF$Z.SCORE = 0
for (cancer in datasets){
  print(cancer)
  genes.not.exprs = c()
  #Read expression data
  norm.counts = read.csv(paste(sep="",output.dir,cancer,"/Normalisation/Norm.Log.Counts.csv"))
  pat.ids.type = gsub("TCGA\\.([[:alnum:]]{2})\\.([[:alnum:]]{4})\\.([[:alnum:]]{2}).*","\\2-\\3", colnames(norm.counts))
  pat.ids = gsub("TCGA\\.([[:alnum:]]{2})\\.([[:alnum:]]{4})\\.([[:alnum:]]{2}).*","\\2", colnames(norm.counts))
  if(cancer=='LGG'){
    pat.ids.norm = paste(pat.ids,'01',sep = '-')
  }else{
  
    pat.ids.norm = paste(pat.ids,'11',sep = '-')
  }

  colnames(norm.counts) = gsub("TCGA\\.([[:alnum:]]{2})\\.([[:alnum:]]{4})\\.([[:alnum:]]{2}).*","\\2-\\3",
                                   colnames(norm.counts))
  cancer.counts = norm.counts[,-which(colnames(norm.counts)%in%pat.ids.norm)]
  colnames(cancer.counts) = gsub("([[:alnum:]]{4})\\-([[:alnum:]]{2}).*","\\1", colnames(cancer.counts))

  RMAF %>% filter(CANCER_TYPE==cancer) -> RMAF.cancer
  for (g in unique(RMAF.cancer$ENTREZ_GENE_ID)){

    if(!(g%in%cancer.counts$X)){
      #print(paste('Check',g,'in',cancer))
      genes.not.exprs = c(g,genes.not.exprs)
      next
    }
    RMAF.cancer %>% filter(ENTREZ_GENE_ID==g) -> aux1

    for (p in unique(aux1$PATIENT_ID)){
      #get all patients same cancer and gender
      if(!(p%in%colnames(cancer.counts))){
        next
      }
      gender = clinical.data[which(clinical.data$PATIENT_ID==p),'GENDER']
      clinical.data %>% filter(CANCER_TYPE==cancer,GENDER==gender,PATIENT_ID%in%colnames(cancer.counts)) -> patients
      exprs = cancer.counts[which(cancer.counts$X == g),patients$PATIENT_ID]
      std = sd(exprs)
      m = mean(as.numeric(exprs))
      z = exprs[p] - m / std
      #Record Z-score
      RMAF[which(RMAF$PATIENT_ID==p&RMAF$ENTREZ_GENE_ID==g),'Z.SCORE'] = z


    }
  }
  RMAF[which(RMAF$CANCER_TYPE==cancer&RMAF$ENTREZ_GENE_ID%in%genes.not.exprs),'EXPRS.GENE'] = 'NOT.EXPRS'
}

###############RESET###################
write.csv(RMAF,paste(output.dir,'AllDataSets/Computed.All.RMAF.csv',sep = ''))

P.Coding.X.genes = read.csv("~/Documents/PhD/GenderAnalysis/Genes.Coding.X.csv")
RMAF = read.csv(paste(output.dir,'AllDataSets/Computed.All.RMAF.csv',sep = ''))
RMAF %>% filter(ENTREZ_GENE_ID%in%P.Coding.X.genes$entrezgene) -> RMAF
RMAF$POLYPHEN_DISCRETE = gsub("(.*)\\(.*\\)","\\1", RMAF$POLYPHEN)
RMAF$SIFT_DISCRETE = gsub("(.*)\\(.*\\)","\\1", RMAF$SIFT)


#Clean low reads
RMAF$SUM = RMAF$REF_COUNT + RMAF$ALT_COUNT

RMAF %>% mutate(RMAF = ifelse(SUM>=10,ALT_COUNT/SUM,0.0)) %>% filter(!(EXPRS.GENE=='NOT.EXPRS'&SUM<10)) %>%
  filter(!(SUM<10&Z.SCORE>-2))-> RMAF

RMAF[which(RMAF$Z.SCORE<=-4&RMAF$SUM<10),'RMAF'] = 1.0
RMAF$LOGIT.RMAF = logit(x = RMAF$RMAF,min = 0-0.1 , max = 1+0.1)
#RMAF[which(RMAF$Z.SCORE<=-4&RMAF$SUM<5),'RMAF'] = 1.0
#DoDensityPlotGenderRMAF(dat = RMAF,dat.col = 'RMAF',fill.col = 'GENDER',title.txt = 'RMAF')
#DoDensityPlotGenderRMAF(dat = RMAF,dat.col = 'Z.SCORE',fill.col = 'GENDER',title.txt = 'Z-Scores')
#Define Expressed Mutation RMAF
#RMAF$EM.RMAF = ifelse(RMAF$RMAF>0.3,1,0)
#RMAF$SM.RMAF = ifelse(RMAF$RMAF<=0.3,1,0)

RMAF$EM.RMAF = ifelse(RMAF$RMAF>=0.75,1,0)
RMAF$SM.RMAF = ifelse(RMAF$RMAF<=0.20,1,0)


#Use genes mutated  in 10+ patients
RMAF %>% group_by(HUGO_SYMBOL, PATIENT_ID) %>% dplyr::summarise(s=max(RMAF)) %>% 
  dplyr::mutate(U_ID = paste(HUGO_SYMBOL, PATIENT_ID,s,sep = '.')) -> aux

RMAF$U_ID = paste(RMAF$HUGO_SYMBOL, RMAF$PATIENT_ID,RMAF$RMAF,sep = '.')
RMAF %>% filter(U_ID %in% aux$U_ID) -> RMAF

# RMAF %>% group_by(HUGO_SYMBOL, PATIENT_ID) %>% dplyr::summarise(s=max(START_POSITION)) %>% 
#   dplyr::mutate(U_ID = paste(HUGO_SYMBOL, PATIENT_ID,s,sep = '.')) -> aux
# 
# RMAF$U_ID = paste(RMAF$HUGO_SYMBOL, RMAF$PATIENT_ID,RMAF$START_POSITION,sep = '.')
# RMAF %>% filter(U_ID %in% aux$U_ID) -> RMAF

RMAF %>% group_by(HUGO_SYMBOL,EM.RMAF,SM.RMAF) %>% dplyr::summarise(EM.SUM=sum(EM.RMAF),SM.SUM=sum(SM.RMAF)) %>% 
  filter(EM.SUM>=10|SM.SUM>=10) -> aux
RMAF %>% filter(HUGO_SYMBOL %in% aux$HUGO_SYMBOL) -> RMAF
rm(aux)

# delt = data.frame()
# for (i in (1:100)){
#   print (i)
#   RMAF %>% mutate(RMAF = ifelse(SUM>=i,ALT_COUNT/SUM,0.0))%>%filter(SUM>=i|(Z.SCORE<=-4))  -> RMAF.aux
#   RMAF.aux[which(RMAF.aux$Z.SCORE<=-4&RMAF.aux$SUM<=i),'RMAF'] = 1.0
#   #RMAF.aux$RMAF = RMAF.aux$ALT_COUNT/RMAF.aux$SUM
#   RMAF.aux$EM.RMAF = ifelse(RMAF.aux$RMAF>=0.5,1,0)
#   RMAF.aux$SM.RMAF = ifelse(RMAF.aux$RMAF<=0.1,1,0)
#   RMAF.aux %>% filter(EM.RMAF==1)%>%group_by(GENDER) %>% dplyr::summarise(n=n()) -> aux.EM
#   RMAF.aux %>% filter(SM.RMAF==1)%>%group_by(GENDER) %>% dplyr::summarise(n=n()) -> aux.SM
#   
#   delt.aux = data.frame(DEPTH = i, FEM.EM = aux.EM[1,2],MAL.EM=aux.EM[2,2], FEM.SM = aux.SM[1,2], MAL.SM=aux.SM[2,2])
#   colnames(delt.aux) =c('DEPTH','FEM.EM','MAL.EM','FEM.SM','MAL.SM')
#   delt = rbind(delt,delt.aux)
# }
# 
# #Check in plot 
# len = 100
# xrange = range(0:len)
# yrange = range(min(delt[2:len,2:5]):max(delt[2:len,2:5]))
# plot(xrange,yrange, type='n')
# colors = rainbow(4)
# linetype = c(1:4)
# 
# for(i in 2:5)
# {
#   lines(delt[2:len,1],delt[2:len,i], type = 'b', lwd=1.5, lty=linetype[i-1],col=colors[i-1])
# }
# legend(6,3200, colnames(delt[,2:5]),col = colors, lty=linetype)

#RMAF %>% filter(SUM>10) -> RMAF
#RMAF$RMAF = RMAF$ALT_COUNT/RMAF$SUM
#RMAF %>% filter(EXPRS.GENE=='EXPRS') -> RMAF



################################################################################################
###############Build SM and EM Matrix for PS tests#############################
################################################################################################
backup = RMAF
RMAF = backup

#RMAF = filter(RMAF,POLYPHEN_DISCRETE=="possibly_damaging"|POLYPHEN_DISCRETE== "probably_damaging")
#RMAF = filter(RMAF,POLYPHEN_DISCRETE == 'benign')
badchars <- "[\xb5]|[\n]|[,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[/]|[\\]|[.]|[_]|[ ]"
clinical.data$PATHOLOGIC_STAGE <- toupper(gsub(pattern=badchars, replacement=".", x=clinical.data$PATHOLOGIC_STAGE))
PS.SM.mat = data.frame(matrix(data = 0,nrow = length(unique(DNA.VAF$PATIENT_ID)), 
                       ncol = length(unique(RMAF$HUGO_SYMBOL))),
                       stringsAsFactors = F, row.names = unique(DNA.VAF$PATIENT_ID))
colnames(PS.SM.mat) = unique(RMAF$HUGO_SYMBOL)
for (g in unique(RMAF$HUGO_SYMBOL)){
  PS.SM.mat[which(rownames(PS.SM.mat)%in%(RMAF[which(RMAF$SM.RMAF==1&RMAF$HUGO_SYMBOL==g),'PATIENT_ID'])),g]  = 1
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
res = DoStatTest(clinical.data = annot,molecular.data = PS.SM.mat,default.confounding.factors = cfd, extra.cofounding.factors = c('SMOKING_STATUS','RACE'))
res$DATASET = 'ALL'
res = as.data.frame(res)
#Add SM numbers and total numbers 
tot.fem = length(rownames(PS.SM.mat)[which(rownames(PS.SM.mat)%in%filter(clinical.data,GENDER=='female')$PATIENT_ID)])
tot.male = length(rownames(PS.SM.mat)[which(rownames(PS.SM.mat)%in%filter(clinical.data,GENDER=='male')$PATIENT_ID)])
res$total.samples = paste(tot.male,'M',":","F",tot.fem)
res$total.males = tot.male
res$total.females = tot.fem
res$silent.mutations = ''
res$silent.mutations.males = 0
res$silent.mutations.females = 0
for(g in colnames(PS.SM.mat)){
  sm.fem = sum(PS.SM.mat[which(rownames(PS.SM.mat)%in%filter(clinical.data,GENDER=='female')$PATIENT_ID),g],na.rm = T)
  sm.male = sum(PS.SM.mat[which(rownames(PS.SM.mat)%in%filter(clinical.data,GENDER=='male')$PATIENT_ID),g],na.rm = T)
  res[which(res$feature==g),'silent.mutations'] = paste(sm.male,'M',":","F",sm.fem)
  res[which(res$feature==g),'silent.mutations.males'] = sm.male
  res[which(res$feature==g),'silent.mutations.females'] = sm.fem
  
}

for (cancer in datasets){
  print(cancer)
  
  ###Filter for Cancer type and genes with at least 10 muts
  RMAF %>% filter(SM.RMAF==1,CANCER_TYPE==cancer) -> RMAF.SM
  RMAF.SM %>% group_by(HUGO_SYMBOL,SM.RMAF) %>% dplyr::summarise(SM.SUM=sum(SM.RMAF)) %>%
    filter(SM.SUM>=10) -> aux
  RMAF.SM %>% filter(HUGO_SYMBOL %in% aux$HUGO_SYMBOL) -> RMAF.SM
  
  DNA.VAF %>% filter(CANCER_TYPE==cancer) -> DNA.cancer
  PS.SM.mat = data.frame(matrix(data = 0,nrow = length(unique(DNA.cancer$PATIENT_ID)), 
                                ncol = length(unique(RMAF.SM$HUGO_SYMBOL))),
                         stringsAsFactors = F, row.names = unique(DNA.cancer$PATIENT_ID))
  colnames(PS.SM.mat) = unique(RMAF.SM$HUGO_SYMBOL)
  for (g in unique(RMAF.SM$HUGO_SYMBOL)){
    RMAF.SM %>% filter(HUGO_SYMBOL==g) -> RMAF.g
    PS.SM.mat[which(rownames(PS.SM.mat)%in%RMAF.g$PATIENT_ID),g] = 1
  }
  clinical.data %>% filter(CANCER_TYPE==cancer) -> clinical.cancer
  annot = clinical.cancer[which(clinical.cancer$PATIENT_ID%in%rownames(PS.SM.mat)),unique(c('PATIENT_ID',cf,cfd))]
  keep = rownames(PS.SM.mat)%in%clinical.cancer$PATIENT_ID
  PS.SM.mat = PS.SM.mat[keep,]
  res.aux = DoStatTest(clinical.data = annot,molecular.data = PS.SM.mat,default.confounding.factors = cfd, extra.cofounding.factors = cf)
  res.aux = as.data.frame(res.aux)
  res.aux$DATASET = cancer
  ####ADD TOTAL NUMBERS
  tot.fem = length(rownames(PS.SM.mat)[which(rownames(PS.SM.mat)%in%filter(clinical.data,GENDER=='female')$PATIENT_ID)])
  tot.male = length(rownames(PS.SM.mat)[which(rownames(PS.SM.mat)%in%filter(clinical.data,GENDER=='male')$PATIENT_ID)])
  res.aux$total.samples = paste(tot.male,'M',":","F",tot.fem)
  res.aux$total.males = tot.male
  res.aux$total.females = tot.fem
  res.aux$silent.mutations.males = 0
  res.aux$silent.mutations.females = 0
  res.aux$silent.mutations = ''
  for(g in colnames(PS.SM.mat)){
    sm.fem = sum(PS.SM.mat[which(rownames(PS.SM.mat)%in%filter(clinical.data,GENDER=='female')$PATIENT_ID),g])
    sm.male = sum(PS.SM.mat[which(rownames(PS.SM.mat)%in%filter(clinical.data,GENDER=='male')$PATIENT_ID),g])
    res.aux[which(res.aux$feature==g),'silent.mutations'] = paste(sm.male,'M',":","F",sm.fem)
    res.aux[which(res.aux$feature==g),'silent.mutations.males'] = sm.male
    res.aux[which(res.aux$feature==g),'silent.mutations.females'] = sm.fem
  }
  
  ###Merge data.frame
  res = rbind(res,as.data.frame(res.aux))
}

write.csv(res,paste(output.dir,'AllDataSets/','SM.Chi.Test.Propensity.KIRC.csv',sep = ''),row.names = F)


###Do p53 interactors tests
res$z.score = sign(res$coef)*sqrt(res$chisq)
res$TP53.INT.LOW = 'No'
res[which(res$feature%in%p53.low),'TP53.INT.LOW'] = 'Yes'
test.tbl = NULL


for (cancer in c('ALL')){
  tmp = filter(res,DATASET==cancer)
  #test.high = t.test(tmp$z.score~tmp$TP53.INT.HIGH)
  #test.med =  t.test(tmp$z.score~tmp$TP53.INT.MED)
  test.low =  wilcox.test(tmp$z.score~tmp$TP53.INT.LOW)
  if (is.null(test.tbl)){
    test.tbl = data.frame(cancer,test.low$statistic,test.low$p.value,'Low.conf')
    colnames(test.tbl) = c("DATASET",'t-statistic','p-value','Interactors')
  }
  else{
    test.tmp = data.frame(cancer,test.low$statistic,test.low$p.value,'Low.conf')
    colnames(test.tmp) = c("DATASET",'t-statistic','p-value','Interactors')
    test.tbl = rbind(test.tbl,test.tmp)
  }
  
  # test.tmp = data.frame(cancer,test.med$statistic,test.med$p.value,'Med.conf')
  # colnames(test.tmp) = c("DATASET",'t-statistic','p-value','Interactors')
  # test.tbl = rbind(test.tbl,test.tmp)
  # test.tmp = data.frame(cancer,test.low$statistic,test.low$p.value,'Low.conf')
  # colnames(test.tmp) = c("DATASET",'t-statistic','p-value','Interactors')
  # test.tbl = rbind(test.tbl,test.tmp)
  
}
write.csv(test.tbl,paste(output.dir,'AllDataSets/','p53.int.SM.wilcox.test.csv',sep = ''),row.names = F)


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
  ###Filter for Cancer type and genes with at least 10 muts
  RMAF %>% filter(EM.RMAF==1,CANCER_TYPE==cancer) -> RMAF.EM
  RMAF.EM %>% group_by(HUGO_SYMBOL,EM.RMAF) %>% dplyr::summarise(EM.SUM=sum(EM.RMAF)) %>% 
    filter(EM.SUM>=1) -> aux
  RMAF.EM %>% filter(HUGO_SYMBOL %in% aux$HUGO_SYMBOL) -> RMAF.EM
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
  PS.EM.mat$PATIENT_ID = rownames(PS.EM.mat)
  PS.EM.mat = PS.EM.mat[keep,]
  PS.EM.mat$PATIENT_ID=NULL
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
    em.fem = sum(PS.EM.mat[which(rownames(PS.EM.mat)%in%filter(clinical.data,GENDER=='female')$PATIENT_ID),g])
    em.male = sum(PS.EM.mat[which(rownames(PS.EM.mat)%in%filter(clinical.data,GENDER=='male')$PATIENT_ID),g])
    res.aux[which(res.aux$feature==g),'expressed.mutations'] = paste(em.male,'M',":","F",em.fem)
    res.aux[which(res.aux$feature==g),'expressed.mutations.females'] = em.fem
    res.aux[which(res.aux$feature==g),'expressed.mutations.males'] = em.male
  }
  
  ###Merge data.frame
  res = rbind(res,as.data.frame(res.aux))
  
}

write.csv(res,paste(output.dir,'AllDataSets/','EM.Chi.Test.Propensity.KIRC.csv',sep = ''),row.names = F)


###Do p53 interactors tests
###Do p53 interactors tests
res$z.score = sign(res$coef)*sqrt(res$chisq)
res$TP53.INT.LOW = 'No'
res[which(res$feature%in%p53.low),'TP53.INT.LOW'] = 'Yes'
test.tbl = NULL



for (cancer in c('ALL',datasets)){
  tmp = filter(res,DATASET==cancer)
  #test.high = t.test(tmp$z.score~tmp$TP53.INT.HIGH)
  #test.med =  t.test(tmp$z.score~tmp$TP53.INT.MED)
  test.low = tryCatch({
    wilcox.test(tmp$z.score~tmp$TP53.INT.LOW)  
  }, error = function(e){NA})
  
  if (is.null(test.tbl)){
    if(is.na(test.low)){
      test.tbl = data.frame(cancer,NA,NA,'Low.conf')
    }else{
      test.tbl = data.frame(cancer,test.low$statistic,test.low$p.value,'Low.conf')  
    }
    
    colnames(test.tbl) = c("DATASET",'t-statistic','p-value','Interactors')
  }
  else{
    if(is.na(test.low)){
      test.tmp = data.frame(cancer,NA,NA,'low.conf')  
    }else{
      test.tmp = data.frame(cancer,test.low$statistic,test.low$p.value,'low.conf')  
    }
    
    colnames(test.tmp) = c("DATASET",'t-statistic','p-value','Interactors')
    test.tbl = rbind(test.tbl,test.tmp)
  }
  
  # test.tmp = data.frame(cancer,test.med$statistic,test.med$p.value,'Med.conf')
  # colnames(test.tmp) = c("DATASET",'t-statistic','p-value','Interactors')
  # test.tbl = rbind(test.tbl,test.tmp)
  # test.tmp = data.frame(cancer,test.low$statistic,test.low$p.value,'Low.conf')
  # colnames(test.tmp) = c("DATASET",'t-statistic','p-value','Interactors')
  # test.tbl = rbind(test.tbl,test.tmp)
  
}
write.csv(test.tbl,paste(output.dir,'AllDataSets/','p53.int.EM.wilcox.test.csv',sep = ''),row.names = F)



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

# Yuan Yuan, Hu Chen
# This script shows you how to perform propensity score weighting and check balance afterwards


library(dplyr)
rm(list=ls())
badchars <- "[\xb5]|[\n]|[,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[/]|[\\]|[.]|[_]|[ ]"
datasets = c("ESCA","HNSC","LUSC","BLCA","LIHC","STAD","LGG","COAD","PAAD","READ","SKCM")
#datasets = "SKCM"
analysis="gender"
scripts.dir="~/Documents/PhD/GenderAnalysis/TCGA/Analysis/PropensityScores/Output" 
clinical.file = "~/Documents/PhD/GenderAnalysis/TCGA/Analysis/all.TCGA.curated.clinical.csv"
out.dir = "~/Documents/PhD/GenderAnalysis/TCGA/Analysis/TCGAExpressionExplorerOutput/"
ps.dir = "/PS/"
ps.file = "PS.csv"
cofounding.factors.ind = c('SMOKING_STATUS','RACE','PATHOLOGIC_STAGE')
cofounding.factors.all = c('SMOKING_STATUS','RACE')
cofounding.factors.default = c('PATIENT_ID','GENDER','AGE_AT_INITIAL_PATHOLOGIC_DIAGNOSIS')
############################
dir.create(file.path(scripts.dir),showWarnings = FALSE)



##### Read clinical data
# headers: source, sample, histologic_type, birth_country, race, gender, age_at_diagnosis, smoking_history, alcohol_history, others
# You may include other factors which you think are confounders
clinical_data <- read.csv(clinical.file, stringsAsFactors=FALSE)
clinical_data_backup = clinical_data
#FILTER NAs and cofounding factors with too many NAs
for (i in c(cofounding.factors.default,cofounding.factors.all))
{
  nr.patients = dim(clinical_data)[1]
  cutoff = ceiling(nr.patients*.25)
  nas = length(clinical_data[which(is.na(clinical_data[,i])),i])
  if(nas>cutoff){
    cofounding.factors.all = cofounding.factors.all[cofounding.factors.all!= i]
  }
  else{
    clinical_data=clinical_data[which(!is.na(clinical_data[,i])),]
    
  }
}

cancer = 'all'
clinical_data = clinical_data[,c(cofounding.factors.default,cofounding.factors.all)]
clinical_data$GENDER <- ifelse(clinical_data$GENDER=="female",1,0)
colnames(clinical_data)[which(colnames(clinical_data)=='GENDER')] <- "Z"
rownames(clinical_data) = clinical_data$PATIENT_ID
clinical_data$PATIENT_ID = NULL

#Calculate Scores for all samples 
# convert clinical variables to dummy
library(dummies)
dummy.feature <- setdiff(colnames(clinical_data),c("Z",'AGE_AT_INITIAL_PATHOLOGIC_DIAGNOSIS'))
clinical_data <- dummy.data.frame(clinical_data, names=cofounding.factors.all)
dummy.list <- attr(clinical_data,"dummies")
rm.col <- c()
for (i in 1:length(dummy.list))
{
  rm.col <- c(rm.col, dummy.list[[i]][length(dummy.list[[i]])])
}
clinical_data <- clinical_data[,-rm.col]
clinical_data$X0 <- rep(1, nrow(clinical_data))
exclude.col <- match(c("Z","X0"), colnames(clinical_data))
colnames(clinical_data) <- toupper(gsub(pattern=badchars, replacement=".", x=colnames(clinical_data)))
#colnames(clinical_data) <- gsub(" ", ".", colnames(clinical_data))
form <- as.formula(paste("Z~",paste(colnames(clinical_data)[-exclude.col],collapse="+"),sep=""))
##### End of processing clinical data

############################
##### perform calculation
source("cal.R")

#Mutation.pri will be the input for the analysis. Each line is a sample while each column is a gene
dir.create(file.path(paste(scripts.dir, "/",cancer,"_",analysis,sep="")),showWarnings = FALSE)
wt.all <- weight.test(clinical_data, form, NULL, is.continuous=TRUE,weight=ifelse(analysis=="gender","MW","ATT"),mirror.plot=TRUE, cancer, "RMAF", outdir=paste(scripts.dir, "/",cancer,"_",analysis,sep=""))
clinical_data$WEIGHTS = wt.all
clinical_data$PATIENT_ID = rownames(clinical_data)
#Save PS for all
dir.create(file.path(paste(out.dir,'AllDataSets',ps.dir,sep = '')),showWarnings = FALSE)
write.csv(clinical_data[,c("PATIENT_ID","WEIGHTS")],paste(out.dir,'AllDataSets',ps.dir,ps.file, sep = ''), row.names = F)


###Calculate for individual Cancers
for (c in datasets)
{
  print(c)
  clinical_data = clinical_data_backup[which(clinical_data_backup$CANCER_TYPE==c),]
  #FILTER NAs and cofounding factors with too many NAs
  cf = cofounding.factors.ind
  for (i in c(cf,cofounding.factors.default))
  {
    nr.patients = dim(clinical_data)[1]
    cutoff = ceiling(nr.patients*.1)
    nas = length(clinical_data[which(is.na(clinical_data[,i])),i])
    if(nas>cutoff){
      cf = cf[cf!= i]
    }
    else{
      clinical_data=clinical_data[which(!is.na(clinical_data[,i])),]
      
    }
  }

  cancer = c
  clinical_data = clinical_data[,c(cofounding.factors.default,cf)]
  clinical_data$GENDER <- ifelse(clinical_data$GENDER=="female",1,0)
  colnames(clinical_data)[which(colnames(clinical_data)=='GENDER')] <- "Z"
  rownames(clinical_data) = clinical_data$PATIENT_ID
  clinical_data$PATIENT_ID = NULL
  
  #Calculate Scores for all samples 
  # convert clinical variables to dummy
  library(dummies)
  dummy.feature <- setdiff(colnames(clinical_data),c("Z",'AGE_AT_INITIAL_PATHOLOGIC_DIAGNOSIS'))
  clinical_data <- dummy.data.frame(clinical_data, names=cf)
  dummy.list <- attr(clinical_data,"dummies")
  rm.col <- c()
  for (i in 1:length(dummy.list))
  {
    rm.col <- c(rm.col, dummy.list[[i]][length(dummy.list[[i]])])
  }
  clinical_data <- clinical_data[,-rm.col]
  clinical_data$X0 <- rep(1, nrow(clinical_data))
  exclude.col <- match(c("Z","X0"), colnames(clinical_data))
  colnames(clinical_data) <- toupper(gsub(pattern=badchars, replacement=".", x=colnames(clinical_data)))
  #colnames(clinical_data) <- gsub(" ", ".", colnames(clinical_data))
  form <- as.formula(paste("Z~",paste(colnames(clinical_data)[-exclude.col],collapse="+"),sep=""))
  ##### End of processing clinical data
  
  ############################
  ##### perform calculation
  source("cal.R")
  
  #Mutation.pri will be the input for the analysis. Each line is a sample while each column is a gene
  dir.create(file.path(paste(scripts.dir, "/",cancer,"_",analysis,sep="")),showWarnings = FALSE)
  wt.all <- weight.test(clinical_data, form, NULL, is.continuous=TRUE,weight=ifelse(analysis=="gender","MW","ATT"),mirror.plot=TRUE, cancer, "RMAF", outdir=paste(scripts.dir, "/",cancer,"_",analysis,sep=""))
  clinical_data$WEIGHTS = wt.all
  clinical_data$PATIENT_ID = rownames(clinical_data)
  #Save PS for cancer
  dir.create(file.path(paste(out.dir,c,ps.dir,sep = '')),showWarnings = FALSE)
  write.csv(clinical_data[,c("PATIENT_ID","WEIGHTS")],paste(out.dir,c,ps.dir,ps.file, sep = ''), row.names = F)
  
    
}



# #rmaf.sum <- summarize.fdr(mutation.pri, rmaf.result, print=TRUE,cutoff=0.05)
# 
# # mRNAseq
# # Read your data. Data format should be the same as mutation.pri
# #mRNAseq.result <- weight.test(clinical_data, form, mRNAseq.pri,is.continuous=TRUE, weight=ifelse(analysis=="gender","MW","ATT"), mirror.plot=TRUE, cancer,"mRNAseq",outdir=paste(scripts.dir, "/",cancer,"_",analysis,sep=""))
# #sum.mRNAseq <- summarize.fdr(mRNAseq.pri, mRNAseq.result)
# 
# # Save the results
# save(mut.result, mRNAseq.result, file=paste(cancer,"_", analysis,"_result.RData",sep=""))
# save(sum.mut, sum.mRNAseq, file=paste(cancer,"_", analysis,"_summary.RData",sep=""))
# 
# ## write summary (significant results) to file
# write.summary(sum.mRNAseq, cancer, analysis,"mRNAseq")
# write.summary(rmaf.sum, cancer, analysis,"rmaf")
# 
# ## write results to file
# write.result(mRNAseq.result, cancer, analysis,"mRNAseq")
# write.result(mut.result, cancer, analysis,"mut")
# 
# 
# ############################
# ##### perform permutation test
# library(doMC)
# library(foreach)
# registerDoMC(15)
# perm.result <- foreach (seed =1:100, .combine='rbind') %dopar%
#     {
#         perm.sum.mut <- c()
#         perm.sum.mRNAseq <- c()
#         
#         ##mRNAseq, which contains continuous value
#         perm.mRNAseq.result <- weight.test(clinical_data, form, mRNAseq.pri,is.continuous=TRUE,  weight=ifelse(analysis=="gender","MW","ATT"), mirror.plot=TRUE, cancer,"mRNAseq",outdir=paste(scripts.dir, "/",cancer,"_",analysis,sep=""), perm=TRUE, seed=seed)
#         perm.sum.mRNAseq <- summarize.fdr(mRNAseq.pri, perm.mRNAseq.result)
# 
#         ## mutation, which contains binary value
#         perm.mut.result <- weight.test(clinical_data, form, mutation.pri, is.continuous=FALSE,weight=ifelse(analysis=="gender","MW","ATT"),mirror.plot=TRUE, cancer, "mut", outdir=paste(scripts.dir, "/",cancer,"_",analysis,sep=""),perm=TRUE, seed=seed)
#         perm.sum.mut <- summarize.fdr(mutation.pri, perm.mut.result)
#         
#         write(c(seed, perm.sum.mRNAseq$n.sig, perm.sum.mut$n.sig), file=paste(scripts.dir, "/",cancer,"_",analysis,"/perm_sig_count_",seed,".txt", sep=""), ncolumns=7)
#         save(seed, perm.mRNAseq.result, perm.mut.result, file=paste(scripts.dir, "/",cancer,"_",analysis,"/perm_sig_result_",seed,".RData", sep=""))
# }
# cutoff <- 0.05
# seedV <- 1:100
# perm.cal(cancer, analysis, "mut", mutation.pri, cutoff=cutoff, seedV=seedV)
# perm.cal(cancer, analysis, "mRNAseq", cnv.pri, cutoff=cutoff, seedV=seedV)
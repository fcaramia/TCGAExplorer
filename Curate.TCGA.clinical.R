library(dplyr)
rm(list=ls())
#Datasets to use
datasets = c("KIRC","LUAD","SARC","ESCA","HNSC","LUSC","BLCA","LIHC","STAD","LGG","COAD","PAAD","READ","SKCM","LUAD",'BRCA','OV','PRAD')
ffpe_tbl = data.table::fread('~/Documents/PhD/Data/TCGA_2016_01_28_BROAD/FFPE.RNA.Files/ffpe_rna_files.csv')
dir = "~/Documents/PhD/Data/TCGA_2016_01_28_BROAD/Clinical.Data/"
data = NULL
badchars <- "[\xb5]|[\n]|[,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[/]|[\\]|[.]|[_]|[ ]"
key_words = c('status','gender','radiation','metas',"age","race","ethnicity","smoke","cigarette","tobacco","stage","histology","alcohol","score",'pam50','pam','ffpe')
key_words = toupper(key_words)
cols = NULL
for (i in datasets)
{
  aux_data = read.delim(paste(dir,i,"/All_CDEs.txt",sep = ""), as.is = T)
  aux_data = t(aux_data)
  colnames(aux_data) = aux_data[1,]
  aux_data = aux_data[-1,]
  rownames(aux_data) <- toupper(rownames(aux_data))
  colnames(aux_data) <- toupper(gsub(pattern=badchars, replacement="_", x=colnames(aux_data)))
  aux_data = as.data.frame(aux_data)
  aux_data$CANCER_TYPE = i
  aux_data$SMOKING_STATUS = NA
  aux_data$PATIENT_ID = toupper(aux_data$PATIENT_ID)
  if ('TOBACCO_SMOKING_HISTORY'%in%colnames(aux_data)){
    aux_data$SMOKING_STATUS = ifelse(aux_data$TOBACCO_SMOKING_HISTORY!='1','yes','no')
  }
   
  if ('RACE'%in%colnames(aux_data)){
    aux_data$RACE = ifelse(aux_data$RACE=='white','white','nonwhite')
  }
 
  for (w in key_words)
  {
    for(c in colnames(aux_data)){
      if(grepl(w,c)){
        cols = unique(c(cols,c))
      }  
    }
    
  }
  if (is.null(data)){
    data = aux_data
  } 
  else{
    cc <- Reduce(intersect,list(colnames(aux_data),colnames(data)))
    data <- bind_rows(data,aux_data)
    data <- data[!duplicated(data),]
    #data <- rbind(data[ , cc, drop=FALSE], aux_data[ , cc, drop=FALSE])    
  }
  
  rm(aux_data)
}


PATIENT_ID = gsub("TCGA\\-([[:alnum:]]{2})\\-([[:alnum:]]{4}).*","\\2",ffpe_tbl$case_id)

data$IS_FFPE = ifelse(data$PATIENT_ID%in%PATIENT_ID,T,F)

write.csv(data,"~/Documents/PhD/GenderAnalysis/TCGA/Analysis/all.TCGA.curated.clinical.csv", row.names = F)

# data[,cols] %>%
#   group_by(GENDER)%>%
#   summarise_each(funs(sum(!is.na(.))))%>% 
#   t -> Info_per_sample
# #data[,cols]%>%summarise_each(funs(sum(!is.na(.))))%>%t -> Info_per_sample
# colnames(Info_per_sample) = Info_per_sample[1,]
# Info_per_sample = as.data.frame(Info_per_sample[-1,])
# Info_per_sample$Total = as.numeric(as.vector(Info_per_sample[,1]))+as.numeric(as.vector(Info_per_sample[,2]))
# Info_per_sample$Clinical.Variable = rownames(Info_per_sample)
# write.csv(Info_per_sample,"~/Documents/PhD/GenderAnalysis/TCGA/Analysis/Clinical.Info.Per.Gender.csv", row.names = F)


##DELETE WHEN DONE###
# muts = read.csv("~/Documents/PhD/GenderAnalysis/TCGA/Analysis/all.reduced.TCGA.curated.mutations.csv")
# muts$PATIENT_ID = gsub("TCGA\\-([[:alnum:]]{2})\\-([[:alnum:]]{4}).*","\\2", muts$SAMPLE_ID)
# muts[,c("PATIENT_ID","CANCER_TYPE")]%>%
#   group_by(PATIENT_ID,CANCER_TYPE)%>%
#   summarise(N_MUTS = n()) %>%
#   left_join(data[,c('PATIENT_ID','PATHOLOGIC_STAGE')]) -> muts_per_sample
# 
# muts_per_sample %>%
#   group_by(PATHOLOGIC_STAGE,CANCER_TYPE) %>%
#   summarise(MUTS_STAGE=median(N_MUTS),N_SAMPLE=n()) -> muts_per_sample_per_stage_per_cancer_type
# 
# write.csv(muts_per_sample_per_stage,"~/Documents/PhD/GenderAnalysis/TCGA/Analysis/muts_per_sample_per_stage.csv", row.names = F)
# 
# write.csv(muts_per_sample_per_stage_per_cancer_type,"~/Documents/PhD/GenderAnalysis/TCGA/Analysis/muts_per_sample_per_stage_per_cancer_type.csv", row.names = F)
# 

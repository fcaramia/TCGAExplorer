rm(list = ls())
library(dplyr)
#TCGA merge mutations curations
#Read Mutations
datasets = c("KIRC","LUAD","ESCA","HNSC","LUSC","BLCA","LIHC","STAD","LGG","COAD","PAAD","READ","SKCM","BRCA","OV")
dir = "~/Documents/PhD/Data/TCGA_2016_01_28_BROAD/Mutation.Data/"
badchars <- "[\xb5]|[\n]|[,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[/]|[\\]|[.]|[_]|[ ]"
mutations = NULL

for (i in datasets)
{
  print(i)
  #aux_data = read.delim(paste(dir,i,"/AllRAWMuts.Oncotator.txt",sep = ""), as.is = T)
  aux_data = read.delim(paste(dir,i,"/AllMuts.Oncotator.txt",sep = ""), as.is = T)
  aux_data$Sample_Id=gsub("TCGA\\-([[:alnum:]]{2})\\-([[:alnum:]]{4})\\-([[:alnum:]]{3}).*","TCGA-\\1-\\2-\\3",aux_data$Tumor_Sample_Barcode)
  aux_data$Cancer.Type = i
  aux_data$DNA_ALT_COUNT = NA
  aux_data$DNA_REF_COUNT = NA
  #Make colnames consistent
  colnames(aux_data) <- toupper(gsub(pattern=badchars, replacement="_", x=colnames(aux_data)))
  aux_data = aux_data[which( aux_data$VARIANT_CLASSIFICATION!="Silent"),]
  
  if("T_ALT_COUNT"%in%colnames(aux_data))
  {
    aux_data$DNA_ALT_COUNT = as.integer(gsub("(.*)\\|(.*)","\\1", aux_data$T_ALT_COUNT))
    aux_data$DNA_REF_COUNT = as.integer(gsub("(.*)\\|(.*)","\\1", aux_data$T_REF_COUNT))
  }
  if("I_TUMOR_VAR_READS"%in%colnames(aux_data))
  {
    aux_data$DNA_ALT_COUNT = as.integer(gsub("(.*)\\|(.*)","\\1", aux_data$I_TUMOR_VAR_READS))
    aux_data$DNA_REF_COUNT = as.integer(gsub("(.*)\\|(.*)","\\1", aux_data$I_TUMOR_REF_READS))
      
  }
  if("I_TVARCOV"%in%colnames(aux_data))
  {
    aux_data$DNA_ALT_COUNT = as.integer(gsub("(.*)\\|(.*)","\\1", aux_data$I_TVARCOV))
    aux_data$DNA_REF_COUNT = as.integer(gsub("(.*)\\|(.*)","\\1", aux_data$I_TTOTCOV)) - 
      as.integer(gsub("(.*)\\|(.*)","\\1", aux_data$I_TVARCOV))
  }
  
  if(is.null(mutations))
  {
    mutations = aux_data 
  }
  else
  {
      #Select common columns 
      cc = intersect(colnames(aux_data),colnames(mutations))
      #Merge datasets
      mutations <- rbind(mutations[ , cc, drop=FALSE], aux_data[ , cc, drop=FALSE])
  }
  #clean up
  rm(aux_data)
  
}

#Discard INDELS and Silent (OPTIONAL)

##Summarise mutations

###Clinical Data####
demo = read.csv("~/Documents/PhD/GenderAnalysis/TCGA/Analysis/all.TCGA.curated.clinical.csv", as.is = T)
demo$PATIENT_ID = toupper(demo$PATIENT_ID)

mutations$PATIENT_ID = gsub("TCGA\\-([[:alnum:]]{2})\\-([[:alnum:]]{4}).*","\\2", mutations$SAMPLE_ID)

mutations[,c("HUGO_SYMBOL","PATIENT_ID")] %>% 
  left_join(demo[c('GENDER','PATIENT_ID')]) %>%
  filter(!is.na(GENDER)) %>% 
  group_by(PATIENT_ID) %>% 
  dplyr::summarise(N.MUTS=n()) %>% filter(N.MUTS<15000) -> validated.patients

mutations %>% 
  left_join(demo[c('GENDER','PATIENT_ID')]) %>%
  filter(!is.na(GENDER)) %>% 
  filter(PATIENT_ID%in%validated.patients$PATIENT_ID) -> mutations 

data = mutations
rm(mutations)
#Save results
#write.csv(data,"~/Documents/PhD/GenderAnalysis/TCGA/Analysis/all.TCGA.curated.mutations.csv", row.names = F)
cols = c("HUGO_SYMBOL","ENTREZ_GENE_ID", "CHROMOSOME","START_POSITION","END_POSITION","VARIANT_TYPE",
         "VARIANT_CLASSIFICATION",
         "SAMPLE_ID","CANCER_TYPE","REFERENCE_ALLELE","TUMOR_SEQ_ALLELE1","TUMOR_SEQ_ALLELE2",
         "TUMOR_SAMPLE_BARCODE","MATCHED_NORM_SAMPLE_BARCODE",'DNA_ALT_COUNT','DNA_REF_COUNT')
write.csv(data[,cols],"~/Documents/PhD/GenderAnalysis/TCGA/Analysis/reduced.all.firehose.TCGA.curated.mutations.csv", row.names = F)









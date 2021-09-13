rm(list = ls())
library(dplyr)
#TCGA merge mutations curations
#Read Mutations
datasets = c("ESCA","HNSC","LUSC","BLCA","LIHC","STAD","LGG","COAD","PAAD","READ","SKCM",'KIRC','LUAD','BRCA','OV','PRAD',"GBM")
#datasets = c('PRAD')
dir = "~/Documents/PhD/Data/TCGA_2016_01_28_BROAD/GDC.Merged.Mutations//"
badchars <- "[\xb5]|[\n]|[,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[/]|[\\]|[.]|[_]|[ ]"
mutations = NULL

files = list.files(dir)

for (f in files)
{
  print(f)
  cancer.type = gsub("TCGA\\.([[:alnum:]]*)\\.([[:alnum:]]*)\\..*DR-7.*","\\1",f)
  if (!(cancer.type %in% datasets)){
    next
  }
  caller = gsub("TCGA\\.([[:alnum:]]*)\\.([[:alnum:]]*)\\..*DR-7.*","\\2",f)
  
  aux_data = read.delim(paste(dir,f,sep = ""), as.is = T, comment.char = '#')
  aux_data$Sample_Id=gsub("TCGA\\-([[:alnum:]]{2})\\-([[:alnum:]]{4})\\-([[:alnum:]]{3}).*","TCGA-\\1-\\2-\\3",aux_data$Tumor_Sample_Barcode)
  aux_data$Cancer.Type = cancer.type
  aux_data$Variant.Caller = caller
  aux_data$DNA_ALT_COUNT = NA
  aux_data$DNA_REF_COUNT = NA
  aux_data$Chromosome = gsub("chr([[:alnum:]]*)","\\1",aux_data$Chromosome)
  
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

mutations %>% 
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

data$MUTATION_ID = paste(data$PATIENT_ID,data$HUGO_SYMBOL,data$CHROMOSOME,data$START_POSITION,
                         data$END_POSITION,data$TUMOR_SEQ_ALLELE2,sep = '.')
data$MUTATION_ID_CALLER = paste(data$MUTATION_ID,data$VARIANT_CALLER, sep = '.')

cols = c("HUGO_SYMBOL","ENTREZ_GENE_ID", "CHROMOSOME","START_POSITION","END_POSITION","VARIANT_TYPE",
       "VARIANT_CLASSIFICATION",'GENDER','PATIENT_ID','VARIANT_CALLER',"MUTATION_ID_CALLER",
       "SAMPLE_ID","CANCER_TYPE","REFERENCE_ALLELE","TUMOR_SEQ_ALLELE1","TUMOR_SEQ_ALLELE2",
       "TUMOR_SAMPLE_BARCODE","MATCHED_NORM_SAMPLE_BARCODE",'DNA_ALT_COUNT','DNA_REF_COUNT',
       "CONSEQUENCE","BIOTYPE","SIFT","POLYPHEN",'IMPACT','FILTER','MUTATION_ID')
data%>% group_by(MUTATION_ID) %>% dplyr::summarise(CALLER = max(VARIANT_CALLER)) -> dat

ids = paste(dat$MUTATION_ID,dat$CALLER, sep = '.')

dat = data[which(data$MUTATION_ID_CALLER%in%ids),cols]



write.csv(dat,"~/Documents/PhD/GenderAnalysis/TCGA/Analysis/reduced.all.raw.TCGA.curated.mutations.csv", row.names = F)



data %>% group_by(GENDER,CANCER_TYPE,PATIENT_ID) %>% dplyr::summarise() %>% 
  group_by(GENDER,CANCER_TYPE) %>% dplyr::summarise(N.P = n()) -> patient.numbers






rm(list = ls())
library(dplyr)
library(ggplot2)
library(gtools)

#TCGA merge mutations curations
#Read Mutations
datasets = c("ESCA","HNSC","LUSC","BLCA","LIHC","STAD","LGG","COAD","PAAD","READ","SKCM","KIRC","LUAD","SARC")
#datasets = c("ESCA","HNSC","LUSC","BLCA","LIHC","STAD","LGG","COAD","READ")
dir = "~/Documents/PhD/Data/TCGA_2016_01_28_BROAD/Mutation.Data/"
dir2 = "~/Documents/PhD/GenderAnalysis/TCGA/Analysis/TCGAExpressionExplorerOutput/"
badchars <- "[\xb5]|[\n]|[,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[/]|[\\]|[.]|[_]|[ ]"
mutations = NULL
genes = c('TP53','KRAS','MYC','PIK3CA')
for (i in datasets)
{
  print(i)
  aux_data = read.delim(paste(dir,i,"/AllMuts.Oncotator.txt",sep = ""), as.is = T)
  aux_data$Sample_Id=gsub("TCGA\\-([[:alnum:]]{2})\\-([[:alnum:]]{4})\\-([[:alnum:]]{3}).*","TCGA-\\1-\\2-\\3",aux_data$Tumor_Sample_Barcode)
  aux_data$Cancer.Type = i
  #Make colnames consistent
  colnames(aux_data) <- toupper(gsub(pattern=badchars, replacement="_", x=colnames(aux_data)))
  
  aux_data=aux_data[,c("HUGO_SYMBOL","ENTREZ_GENE_ID","CHROMOSOME","START_POSITION","END_POSITION","VARIANT_TYPE","VARIANT_CLASSIFICATION","SAMPLE_ID","CANCER_TYPE","REFERENCE_ALLELE","TUMOR_SEQ_ALLELE1","TUMOR_SEQ_ALLELE2","TUMOR_SAMPLE_BARCODE","MATCHED_NORM_SAMPLE_BARCODE")]
  
  
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

###Clinical Data####
demo = read.csv("~/Documents/PhD/GenderAnalysis/TCGA/Analysis/all.TCGA.curated.clinical.csv", as.is = T)
demo$PATIENT_ID = toupper(demo$PATIENT_ID)

mutations$PATIENT_ID = gsub("TCGA\\-([[:alnum:]]{2})\\-([[:alnum:]]{4}).*","\\2", mutations$SAMPLE_ID)

##Summarise mutations
mutations %>% 
  left_join(demo[c('GENDER','PATIENT_ID')]) %>%
  filter(!is.na(GENDER)) %>% 
  group_by(PATIENT_ID) %>% 
  dplyr::summarise(N.MUTS=n()) %>% filter(N.MUTS<10000) -> validated.patients
 
mutations %>% 
  left_join(demo[c('GENDER','PATIENT_ID')]) %>%
  filter(PATIENT_ID%in%validated.patients$PATIENT_ID) -> muts 
  

muts %>% filter(VARIANT_CLASSIFICATION!='Silent') %>%
group_by(PATIENT_ID,GENDER,CANCER_TYPE,HUGO_SYMBOL) %>%
  dplyr::summarise(N.MUTS=n()) %>% 
  group_by(GENDER,CANCER_TYPE,HUGO_SYMBOL) %>% 
  dplyr::summarise(N.MUTS=n()) -> all.muts
  
  

samples=NULL
for (cancer in datasets)
{
  norm.counts = read.csv(paste(dir2,cancer,'/Normalisation/Norm.Log.Counts.csv',sep = ""), as.is = T)
  rownames(norm.counts) = norm.counts[,1]
  norm.counts = norm.counts[,-1]
  #########################  
  samp.names = data.frame(SAMPLE_ID=colnames(norm.counts),SAMPLE_CODE = NA, SAMPLE_TYPE = NA,stringsAsFactors=F, 
                             CANCER_TYPE = cancer)
  samp.names$SAMPLE_CODE = gsub("TCGA\\.([[:alnum:]]{2})\\.([[:alnum:]]{4}\\.([[:alnum:]]{2})).*","\\3", samp.names$SAMPLE_ID)
  samp.names$SAMPLE_TYPE <- ifelse(as.integer(samp.names$SAMPLE_CODE)>=20, 'Control',
                              ifelse(as.integer(samp.names$SAMPLE_CODE)>=10, 'Normal', 'Tumor'))
  
  if(is.null(samples))
  {
    samples = samp.names
  }
  else
  {
    samples = rbind(samples,samp.names)
  }
  
  rm(norm.counts,samp.names)
}

samples$PATIENT_ID = gsub("TCGA\\.([[:alnum:]]{2})\\.([[:alnum:]]{4}).*","\\2", samples$SAMPLE_ID)
samples %>% 
  left_join(demo[c('GENDER','PATIENT_ID')]) %>%
  filter(!is.na(GENDER)) -> samples.gender

samples.gender %>%
  group_by(CANCER_TYPE,GENDER,SAMPLE_TYPE) %>%
  dplyr::summarise(N.PATS=n())  -> samp_res


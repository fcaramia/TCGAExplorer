#TCGA merge mutations curations
#Read Mutations
rm(list = ls())
datasets = c("ESCA","HNSC","LUSC","BLCA","LIHC","PAAD","READ","SKCM",'LGG',
             'COAD',"STAD",'LUAD','KIRC.FEMALES','KIRC.MALES','BRCA','OV','PRAD')

#datasets = c("BLCA","LUSC","LUAD","SKCM",'COAD','HNSC',"LGG")

dir="~/Documents/PhD/Data/TCGA_2016_01_28_BROAD/Mutation.Data/"
badchars <- "[\xb5]|[\n]|[,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[/]|[\\]|[.]|[_]|[ ]"

exmut = NULL
for (i in datasets)
{
  print(i)
  aux_data = read.csv(paste(dir,gsub("(.*)\\..*","\\1",i),"/ExMut.",i,".csv",sep = ""), as.is = T)
  aux_data = unique(aux_data)
  
  #Make colnames consistent
  colnames(aux_data) <- toupper(gsub(pattern=badchars, replacement="_", x=colnames(aux_data)))
  if(is.null(exmut))
  {
    exmut = aux_data 
  }
  else
  {
    #Select common columns 
    cc = intersect(colnames(aux_data),colnames(exmut))
    #Merge datasets
    exmut <- rbind(exmut[ , cc, drop=FALSE], aux_data[ , cc, drop=FALSE])
  }
  #clean up
  rm(aux_data)
  
}

#FILTER MUTATIONS to SNP and non Silent
#Mutations already filtered
# mutations = read.csv("~/Documents/PhD/GenderAnalysis/TCGA/Analysis/reduced.all.TCGA.curated.mutations.csv", as.is = T)
# mutations$PATIENT = gsub("TCGA\\-([[:alnum:]]{2})\\-([[:alnum:]]{4}).*","\\2", mutations$TUMOR_SAMPLE_BARCODE)
# mutations$MUTATION_ID = paste(sep = ".",mutations$PATIENT,mutations$CHROMOSOME,mutations$START_POSITION,
#                               mutations$TUMOR_SEQ_ALLELE2)
# exmut$PATIENT = gsub("TCGA\\-([[:alnum:]]{2})\\-([[:alnum:]]{4}).*","\\2", exmut$SAMPLE)
# exmut$MUTATION_ID = paste(sep = ".",exmut$PATIENT,exmut$CHROMOSOME,exmut$START_POSITION,
#                               exmut$TUMOR_SEQ_ALLELE2)

#exmut = exmut[which(exmut$MUTATION_ID%in%mutations$MUTATION_ID),]

##Remove duplicates, keep highest sum
exmut$SUM = exmut$REF_COUNT + exmut$ALT_COUNT

exmut <- exmut[order(exmut$MUTATION_ID, -abs(exmut$SUM) ), ]
exmut = exmut[ !duplicated(exmut$MUTATION_ID), ]
exmut <- exmut[order(exmut$CANCER_TYPE, exmut$PATIENT_ID ), ]

#Save results
write.csv(exmut,"~/Documents/PhD/GenderAnalysis/TCGA/Analysis/all.TCGA.ExMut.csv", row.names = F)
#write.csv(exmut,"~/Desktop/TCGA.ExMut.csv", row.names = F)
#clean up
rm(mutations,exmut)

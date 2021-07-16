library(sevenbridges)
library(data.table)

rm(list = ls())


#######################
germs = fread('~/Documents/PhD/GenderAnalysis/TCGA/Analysis/Germline.Normals.ExMut.csv')
muts.file = 'Germline.Tumors.For.ExMutv5.AllChromosomes.csv'

somat = fread('~/Documents/PhD/GenderAnalysis/TCGA/Analysis/Mutations.For.ExMutv5.AllChromosomes.csv')

aux = as.vector(somat[match(germs$SOMATIC,somat$MUTATION_ID),'SAMPLE_ID'])
germs$SOMATIC_SAMPLE_ID = aux$SAMPLE_ID
rm(aux)

germline.samples = unique(germs$SOMATIC_SAMPLE_ID)
###########################


####Create list of samples ids


auth <- Auth(token = "08f0557180254839968586d622b4dfac", platform = "cgc")

project <- auth$project(id = "PAPENFUSS/xchrmexplorer/")

mutfile <- project$file(name = muts.file, exact = TRUE)
#Get a files metadata
#f = project$file(id = "5c38343fe4b08832b6df151f")$meta()
#$sample_id use filter files 


##Get all files
#fs = project$file(complete = TRUE,)

##Get files by metadata
files_ids = list()
for (s in germline.samples){
  print(s)
  Sys.sleep(5)
  small_list = project$file(metadata=list(sample_id=s))
  if(length(small_list)>1){
    for (m in seq(1,length(small_list)))
    {
      name = small_list[[m]]$name
      if (endsWith(name,'bam')){
        files_ids = append(files_ids, small_list[[m]]$id)
      }
    }  
    
  }
  if(length(small_list)==1){
    name = small_list$name
    if (endsWith(name,'bam')){
      files_ids = append(files_ids, small_list$id)
    }
  }
}

bams_list = NULL
for (f in files_ids){
  bam = project$file(id = f)
  bams_list = append(bams_list,bam)
  Sys.sleep(1)
}

##get all tasks
tasks = project$task(complete=TRUE)

##get task by id
example_task = project$task(id = 'c53c6bd3-5b8d-4b0d-be8b-588f36493a97')

inputs = list(
  bam_in = bams_list,
  mutations = mutfile
  )

project$task_add(name = "ExmutWF-v.0.2 run - Germline Tumors All Chrom", app = "PAPENFUSS/xchrmexplorer/exmutwf-v-0-2/32", 
                 use_interruptible_instances = TRUE, inputs = inputs )


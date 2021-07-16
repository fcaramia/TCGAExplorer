rm(list = ls())
library(dplyr)
library(data.table)
library(tidyr)
#TCGA merge mutations curations
#Read Mutations
#datasets = c("KIRC","LUAD","ESCA","HNSC","LUSC","BLCA","LIHC","STAD","LGG","COAD","PAAD","READ","SKCM","BRCA","OV",'PRAD')
datasets = c("LIHC","STAD","LGG","COAD","PAAD","READ","SKCM","BRCA","OV",'PRAD')
dir = "~/Documents/PhD/Data/TCGA_2016_01_28_BROAD/Methylation.Data/"
output_dir = "~/Documents/PhD/GenderAnalysis/TCGA/Analysis/TCGAExpressionExplorerOutput/"

for (i in datasets)
{
  print(i)
  #aux_data = read.delim(paste(dir,i,"/AllRAWMuts.Oncotator.txt",sep = ""), as.is = T)
  tar_file = paste0("gdac.broadinstitute.org_",i,".Merge_methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.Level_3.2016012800.0.0.tar.gz")
  target_file = paste0("gdac.broadinstitute.org_",i,".Merge_methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.Level_3.2016012800.0.0/",i,".methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt")
  print('Decompressing File')
  untar(tarfile = paste0(dir,i,'/20160128/',tar_file),files=target_file, exdir = "/tmp")
  print("Filtering meta data out of file")
  system( paste0("awk \'BEGIN{FS=\"\\t\"};{printf \"%s\\t%s\\t%s\\t%s\\t%s\",$1,$3,$4,$5,$2; for (i=6;i<NF;i+=4) printf \"\\t%s\", $i; printf \"\\n\"}\' ","/tmp/",target_file,"| sed '2d' >/tmp/out_num.txt"))
  
  print("Opening Methyl file")
  #methyl_data_tbl = read.table(file = '/tmp/out_num.txt', sep='\t',header = T, row.names = 1)
  methyl_data_tbl = fread(file = '/tmp/out_num.txt',sep = '\t',header = T, stringsAsFactors = F, data.table = T)
  
  #annot = methyl_data_tbl[,c(1,3,4,5)]
  #colnames(annot) = unlist(annot[1,])
  #annot = annot[-1,]
  print('Methyl file filtering')
  colnames(methyl_data_tbl)[c(1,2,3,4)] = c("Probe", "Gene_Symbol","Chromosome","Coordinate")
  methyl_data_tbl = methyl_data_tbl[complete.cases(methyl_data_tbl), ]
  methyl_data_tbl = separate_rows(methyl_data_tbl,Gene_Symbol,sep = ";")
  
  # Open the gene expression
  # Use raw counts and avoid filtering 
  print("Getting Gene Expression")
  expr_tbl = fread(paste0(output_dir,i,"/Normalisation","/Norm.Log.Counts.csv"))
  # Make gene ids unique 
  expr_tbl$V1 = gsub(".*\\|(.*)","\\1", expr_tbl$V1)
  expr_tbl$sum = rowSums(expr_tbl[,2:dim(expr_tbl)[2]])
  expr_tbl <- expr_tbl[order(expr_tbl$V1, -abs(expr_tbl$sum) ), ]
  expr_tbl = expr_tbl[ !duplicated(expr_tbl$V1), ]
  colnames(expr_tbl)[1] = 'Gene_Symbol'
  
  # Match sample ids in expr and methylation
  colnames(expr_tbl) = gsub("TCGA\\.([[:alnum:]]{2})\\.([[:alnum:]]{4})\\.([[:alnum:]]{2}).*","TCGA-\\2-\\3",
                               colnames(expr_tbl))
  colnames(methyl_data_tbl) = gsub("TCGA-([[:alnum:]]{2})-([[:alnum:]]{4})-([[:alnum:]]{2}).*","TCGA-\\2-\\3",
                            colnames(methyl_data_tbl))
  
  
  # Keep samples in both tables
  sample_lst = intersect(colnames(expr_tbl),colnames(methyl_data_tbl))
  expr_tbl = expr_tbl %>% select(any_of(sample_lst))
  methyl_data_tbl = methyl_data_tbl %>% select(any_of(c("Probe","Chromosome","Coordinate",sample_lst)))
  
  # Keep probes with expressed genes
  methyl_data_tbl = methyl_data_tbl %>% filter(Gene_Symbol%in%expr_tbl$Gene_Symbol)
  methyl_data_tbl = methyl_data_tbl[order(methyl_data_tbl$Gene_Symbol),]
  # Calculate correlations to expression
  print("Calcutating correlations")
  cor_tbl <- apply(X = methyl_data_tbl, MARGIN = 1,
                     FUN = function (x) cor(as.numeric(x[-c(1,2,3,4)]),
                                            as.numeric(expr_tbl[which(expr_tbl[,1]==as.character(x[4])),-1])))

  methyl_data_tbl$expr_cor = cor_tbl
  min_probs_tbl  = methyl_data_tbl %>% 
    group_by(Gene_Symbol) %>% slice(which.min(expr_cor) )
  
  print("Saving")
  head_cols_lst = c("Probe","Gene_Symbol", "Chromosome","Coordinate","expr_cor")
  res = min_probs_tbl[,c(head_cols_lst,colnames(min_probs_tbl)[-match(head_cols_lst,colnames(min_probs_tbl))])]
  dir.create(paste0(output_dir,i,'/Methylation'),showWarnings = F, recursive = T)
  fwrite(res, paste0(output_dir,i,'/Methylation/methylation_estimates.csv'))
  
  #clean up
  print("Cleaning")
  rm(cor_tbl)
  rm(expr_tbl)
  rm(methyl_data_tbl,min_probs_tbl,res, head_cols_lst, sample_lst)
  file.remove(paste0('/tmp/',target_file))
  file.remove(paste0('/tmp/out_num.txt'))
  gc(verbose = T, full = T)
}


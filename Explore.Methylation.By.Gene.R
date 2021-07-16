rm(list = ls())
gc(verbose = T, full = T)
library(dplyr)
library(data.table)
library(tidyr)
library("IlluminaHumanMethylation450kanno.ilmn12.hg19")
source("~/Documents/Work/repositories/loi-lab/code/RNA.Analysis/GeneSetAnalysis.R")
#TCGA merge mutations curations
#Read Mutations
datasets = c("KIRC","LUAD","ESCA","HNSC","LUSC","BLCA","LIHC","STAD","LGG","COAD","PAAD","READ","SKCM","BRCA","OV",'PRAD')
datasets = c("BRCA")
data("IlluminaHumanMethylation450kanno.ilmn12.hg19")
annotation.table = as.data.frame(getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19))
annotation.table = annotation.table %>% 
  dplyr::select(chr,pos,strand,Name,Regulatory_Feature_Group,UCSC_RefGene_Name,UCSC_RefGene_Group) %>%
  mutate(probe = Name, 
         Name= NULL, 
         symbol = UCSC_RefGene_Name, 
         UCSC_RefGene_Name=NULL,
         probe_loc = UCSC_RefGene_Group ,
         UCSC_RefGene_Group = NULL,
         probe_type = Regulatory_Feature_Group,
         Regulatory_Feature_Group = NULL,  
         chr_ensembl = gsub("chr","",chr)) %>%
  #filter(probe_type %in% c("Promoter_Associated","Promoter_Associated_Cell_type_specific")) %>%
  separate_rows(symbol,sep = ";") %>%
  separate_rows(probe_loc,sep = ";") %>%
  distinct() %>%
  group_by_at(vars(-probe_loc)) %>%
  summarise(probe_loc = paste(probe_loc,collapse = ';'))
  
dir = "~/Documents/PhD/Data/TCGA_2016_01_28_BROAD/Methylation.Data/"
promoters_tbl = GetHumanHg19GenePromoter()
target_range_obj <- GRanges(promoters_tbl$chromosome_name,
                            IRanges(start=promoters_tbl$chromosome_start,
                                    end=promoters_tbl$chromosome_end),
                            ensembl_promoter_id = promoters_tbl$regulatory_stable_id)

ranges_obj <- GRanges(annotation.table$chr_ensembl,
                      IRanges(annotation.table$pos, width=1),
                      Probe = annotation.table$probe)
overlap_tbl = as.data.frame(mergeByOverlaps(ranges_obj,target_range_obj))
overlap_tbl = dplyr::select(overlap_tbl,
Probe,
ensembl_promoter_id,
target_range_obj.start,
target_range_obj.end)
colnames(overlap_tbl) = c("probe","promoter_id","promoter_start","promoter_end")
annotation.table = left_join(annotation.table,overlap_tbl) %>% distinct()

## Overlap with TSS
tss_tbl = GetHumanHg19GeneTSS()
target_range_obj <- GRanges(tss_tbl$chromosome_name,
                            IRanges(start=tss_tbl$transcription_start_site-1500,
                                    end=tss_tbl$transcription_start_site+1500),
                            tss =tss_tbl$transcription_start_site)

overlap_tbl = as.data.frame(mergeByOverlaps(ranges_obj,target_range_obj))
overlap_tbl = dplyr::select(overlap_tbl,
                            Probe,
                            tss
                            )
colnames(overlap_tbl) = c("probe","tss")
annotation.table = left_join(annotation.table,overlap_tbl) %>% 
  distinct() %>%
  group_by_at(vars(-tss)) %>%
  summarise(tss_dist = min(abs(tss - pos)))

output_dir = "~/Documents/PhD/GenderAnalysis/TCGA/Analysis/TCGAExpressionExplorerOutput/"
genes =  c("RAB25","TP53","CFH","AHNAK","LAMA4","FGFR3","CD163","COL11A1","COL6A3","IGHG1","CTNNB1","TNC","CEP250","RBM19","MYH9","FN1","COL5A1","PHF3") 
hyper_methyl_brca = c("APC","BRCA1","CDKN2A","GSTP1","CCND2","PTEN","RARB","RASSF1","ZMYND10")

## Read X study genes
x_genes_lst = read.csv("~/Documents/PhD/GenderAnalysis/TCGA/Analysis/TCGAExpressionExplorerOutput/AllDataSets/SM.Chi.Test.Propensity.All.csv")
x_genes_lst = x_genes_lst[which(x_genes_lst$fdr<0.06),'feature']
genes = c(genes,x_genes_lst, hyper_methyl_brca)
chromosomes = c("chrX")
for (i in datasets)
{
  print(i)
  #aux_data = read.delim(paste(dir,i,"/AllRAWMuts.Oncotator.txt",sep = ""), as.is = T)
  tar_file = paste0("gdac.broadinstitute.org_",i,".Merge_methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.Level_3.2016012800.0.0.tar.gz")
  target_file = paste0("gdac.broadinstitute.org_",i,".Merge_methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.Level_3.2016012800.0.0/",i,".methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt")
  print('Decompressing File')
  start_time <- Sys.time()
  untar(tarfile = paste0(dir,i,'/20160128/',tar_file),files=target_file, exdir = "/tmp")
  end_time <- Sys.time()
  print(end_time - start_time)
  start_time <- Sys.time()
  print("Filtering meta data out of file")
  system( paste0("awk \'BEGIN{FS=\"\\t\"};{printf \"%s\\t%s\",$1,$2; for (i=6;i<NF;i+=4) printf \"\\t%s\", $i; printf \"\\n\"}\' ","/tmp/",target_file,"| sed '2d' >/tmp/out_num.txt"))
  # system( paste0("awk \'BEGIN{FS=\"\\t\"};{printf \"%s\\t%s\\t%s\\t%s\",$1,$3,$4,$5; printf \"\\n\"}\' ","/tmp/",target_file,"| sed '2d' >/tmp/out_annot.txt"))
  end_time <- Sys.time()
  print(end_time - start_time)
  print("Opening Methyl file")
  start_time <- Sys.time()
  methyl_data_tbl = fread(file = '/tmp/out_num.txt',sep = '\t',header = T, stringsAsFactors = F, data.table = T)
  end_time <- Sys.time()
  print(end_time - start_time)
  print('Methyl file filtering')
  colnames(methyl_data_tbl)[c(1)] = "probe"
  
  methyl_data_tbl = methyl_data_tbl[complete.cases(methyl_data_tbl), ]
  methyl_data_tbl = methyl_data_tbl %>% left_join(annotation.table)
  
  # Filter genes of interest
  methyl_data_tbl = methyl_data_tbl%>% filter(symbol%in%genes|chr%in%chromosomes)
  
  print("Getting Gene Expression")
  expr_tbl = fread(paste0(output_dir,i,"/Normalisation","/Norm.Log.Counts.csv"))
  # Make gene ids unique 
  expr_tbl$V1 = gsub(".*\\|(.*)","\\1", expr_tbl$V1)
  expr_tbl$sum = rowSums(expr_tbl[,2:dim(expr_tbl)[2]])
  expr_tbl <- expr_tbl[order(expr_tbl$V1, -abs(expr_tbl$sum) ), ]
  expr_tbl = expr_tbl[ !duplicated(expr_tbl$V1), ]
  colnames(expr_tbl)[1] = 'symbol'
  
  # Match sample ids in expr and methylation
  colnames(expr_tbl) = gsub("TCGA\\.([[:alnum:]]{2})\\.([[:alnum:]]{4})\\.([[:alnum:]]{2}).*","TCGA-\\2-\\3",
                            colnames(expr_tbl))
  colnames(methyl_data_tbl) = gsub("TCGA-([[:alnum:]]{2})-([[:alnum:]]{4})-([[:alnum:]]{2}).*","TCGA-\\2-\\3",
                                   colnames(methyl_data_tbl))
  
  
  # Keep samples in both tables
  sample_lst = intersect(colnames(expr_tbl),colnames(methyl_data_tbl))
  expr_tbl = expr_tbl %>% dplyr::select(any_of(sample_lst))
  methyl_data_tbl = methyl_data_tbl %>% 
    dplyr::select(any_of(c("probe",
                    "chr",
                    "pos",
                    "strand",
                    "probe_type",
                    "probe_loc",
                    "promoter_id",
                    "promoter_start",
                    "promoter_end",
                    "tss_dist",
                    sample_lst)))
  
  # Keep probes with expressed genes
  methyl_data_tbl = methyl_data_tbl[order(methyl_data_tbl$symbol),]
  
  # Gather methyl data
  methyl_data_tbl = methyl_data_tbl %>%
    gather(key = Sample_ID, 
           value = Methyl_Beta_Val,
           -probe,
           -symbol,
           -chr, 
           -pos,
           -probe_type,
           -probe_loc,
           -promoter_id,
           -promoter_start,
           -promoter_end,
           -tss_dist,
           -strand)%>%
    mutate(merge_id = paste0(Sample_ID,'-',symbol)) 
  
  # Gather Expression Data
  expr_tbl = expr_tbl %>% 
    filter(symbol%in%methyl_data_tbl$symbol) %>%
    gather(key = Sample_ID, value = Expr_Val, -symbol ) %>%
    mutate(merge_id = paste0(Sample_ID,'-',symbol)) 
  
  #Merge
  methyl_data_tbl = methyl_data_tbl %>%
    left_join(expr_tbl[,c("Expr_Val","merge_id")],by = c('merge_id'="merge_id"))
  
  methyl_data_tbl$merge_id = NULL
  methyl_data_tbl$Expr_Val = ifelse(is.na(methyl_data_tbl$Expr_Val),0,methyl_data_tbl$Expr_Val)
  fwrite(methyl_data_tbl, paste0(output_dir,i,'/Methylation/methylation_expression_gene_set.csv'))
  
  print("Cleaning")
  
  rm(methyl_data_tbl, expr_tbl)
  file.remove(paste0('/tmp/',target_file))
  file.remove(paste0('/tmp/out_num.txt'))
  gc(verbose = T, full = T)

}


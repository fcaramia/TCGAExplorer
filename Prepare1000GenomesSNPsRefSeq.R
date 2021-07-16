rm(list = ls())
library(data.table)
library(tidyr)
library(dplyr)
library(GenomicRanges)
# germs = fread('~/Documents/PhD/GenderAnalysis/TCGA/Analysis/Filtered1000GenomesXchrForExMutv6.csv')
# 
# 
# gtfhg19 = fread("~/Documents/PhD/Data/GTF.hg19/Homo_sapiens.GRCh37.87.gtf",skip = 5)
# germs$merge_id = paste(germs$chr,germs$pos,germs$alt, sep = "_")
# # filter the X
# gtfhg19 = gtfhg19 %>% filter(V1 == "X",V3=="transcript") %>% select(V1,V4,V5) %>% distinct()
# target_range_obj <- GRanges(gtfhg19$V1,
#                             IRanges(start=gtfhg19$V4,
#                                     end=gtfhg19$V5)
#                             )
# 
# ranges_obj <- GRanges(germs$chr,
#                       IRanges(germs$pos, width=1),merge_id = germs$merge_id)
# 
# overlap_tbl = as.data.frame(mergeByOverlaps(ranges_obj,target_range_obj))
# 
# germs =  germs %>% filter(merge_id%in%overlap_tbl$merge_id)
# 
# germs = germs%>% mutate(AF = as.double(gsub(".*;AF=(.*);AN.*","\\1",info))) %>% filter(AF >= 0.2 , AF <=.8 )
# 
# germs = germs %>% select(chr,pos,ref,alt,var,sample_id,AF)
# 
# fwrite(germs,"../Filtered1000GenomesXchrForExMutv6.csv")

## Filter Seven Bridges 10 sample file

rm(list = ls())
germs = fread('~/Documents/PhD/GenderAnalysis/TCGA/Analysis/Filtered1000GenomesXchrSB10Samples.csv')

germs = germs %>% mutate(sum = if_else(ref_alleles>=alt_alleles&alt_alleles>=other_alleles,ref_alleles+alt_alleles,
                                       if_else(ref_alleles>=alt_alleles&other_alleles>=alt_alleles,ref_alleles+other_alleles,
                                               if_else(alt_alleles>=ref_alleles&other_alleles>=ref_alleles,alt_alleles+other_alleles,alt_alleles+ref_alleles))))
#Filter
germs = germs %>% 
  filter(sum > 10) %>% 
  mutate(AF = as.double(gsub(".*;AF=(.*);AN.*","\\1",info))) %>% 
  mutate(rmaf = if_else(ref_alleles>=alt_alleles&alt_alleles>=other_alleles,alt_alleles/sum,
                       if_else(ref_alleles>=alt_alleles&other_alleles>=alt_alleles,other_alleles/sum,
                               alt_alleles/sum))) %>%
  filter(AF>0.05)


# Summarise
germs_sum = germs %>% 
  group_by(pos) %>% 
  summarise(n=n(),rmaf = mean(rmaf)) %>% 
  filter(n>3,rmaf>0.05,rmaf<0.995)

  
germs = germs %>% filter(pos %in% germs_sum$pos) %>%
  select(chr,pos,ref,alt,var,AF) %>% mutate(sample_id = "1000GenomesGold") %>% distinct()


fwrite(germs,"../Filtered1000GenomesXchrForExMutv610SamplesGold.csv")






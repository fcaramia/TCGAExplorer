library(dplyr)
rm(list = ls())
chrom.sizes = read.csv("~/Documents/PhD/GenderAnalysis/TCGA/Analysis/hg19.chrom.sizes.txt",sep = '\t', 
                       col.names = c('CHROM', "LENGTH"))

genes.entrez = read.csv("~/Documents/PhD/Data/gene_db/entrez_hgnc.csv")

gene.lengths = read.csv("~/Documents/PhD/GenderAnalysis/TCGA/Analysis/gene.sizes.csv")

valid.chrm = c('1','2','3','4','5','6','7','8','9','10','11','12','13',
               '14','15','16','17','18','19','20','21','22','X','Y')

gene.lengths %>% 
  filter(Chromosome.scaffold.name%in%valid.chrm) %>%
  filter(!is.na(HGNC.symbol)&HGNC.symbol!='') %>%
  mutate(GENE.LENGTH = Transcript.length..including.UTRs.and.CDS.) %>%
  left_join(genes.entrez,by = c("HGNC.symbol" = "hgnc_symbol")) %>%
  filter(!is.na(entrezgene)) %>%
  group_by(HGNC.symbol,Chromosome.scaffold.name,entrezgene,Gene.type) %>%
  summarise(GENE.LENGTH = max(GENE.LENGTH)) -> res

colnames(res) = c("GENE.SYMBOL","CHROMOSOME","ENTREZ.ID","GENE.TYPE","GENE.LENGTH")

write.csv(res,"~/Documents/PhD/GenderAnalysis/TCGA/Analysis/All.genes.lengths.csv",row.names = F )
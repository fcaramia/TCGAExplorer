rm(list = ls())
library(data.table)

gtf = fread(input = '~/Documents/PhD/Data/GTF.hg19/Homo_sapiens.GRCh37.87.gtf',header = F,skip = '#!')
gtf = gtf[,1:5]
colnames(gtf) = c("CHROMOSOME","SOURCE","TYPE","START","END")

gtf = filter(gtf, CHROMOSOME=='X',TYPE == 'exon')

gtf =  gtf[order(gtf$START),]



library(dplyr)
library(ggplot2)
gtf %>% 
  dplyr::arrange(START) %>% 
  dplyr::group_by(g = cumsum(cummax(dplyr::lag(END, default = dplyr::first(END))) < START)) %>% 
  dplyr::summarise(START = dplyr::first(START), END = max(END)) -> gtf.sorted

gtf$color = "red"

p <- ggplot() 
p + geom_segment(data = gtf,aes(x=START,y=0, xend = END, yend=0,color=color), size=4, show.legend = F)



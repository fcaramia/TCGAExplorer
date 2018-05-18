rm(list = ls())
library(dplyr)

#Read all interactors

all.int = read.csv("~/Documents/PhD/GenderAnalysis/TP53_interactions/String.p53.interactions.csv")

all.int %>% filter(chromosome_name.x == 'X' | chromosome_name.y == 'X') -> all.x.int
all.int %>% filter(chromosome_name.x != 'X' & chromosome_name.y != 'X' &
                   chromosome_name.x != 'Y' & chromosome_name.y != 'Y') -> all.auto.int

all.x.int %>% filter(score>=800) %>% select(hgnc_symbol.x,hgnc_symbol.y) -> x.800
all.x.int %>% filter(score>=500) %>% select(hgnc_symbol.x,hgnc_symbol.y) -> x.500


l.all = c(as.vector(all.x.int$hgnc_symbol.x),as.vector(all.x.int$hgnc_symbol.y))
l.all = paste(unique(l.all[which(l.all!='TP53')]),collapse = ',')

l.all.auto = c(as.vector(all.auto.int$hgnc_symbol.x),as.vector(all.auto.int$hgnc_symbol.y))
l.all.auto.l = paste(unique(l.all.auto[which(l.all.auto!='TP53')]),collapse = ',')
l.all.auto = unique(l.all.auto[which(l.all.auto!='TP53')])



l.800 = c(as.vector(x.800$hgnc_symbol.x),as.vector(x.800$hgnc_symbol.y))
l.800 = paste(unique(l.800[which(l.800!='TP53')]),collapse = ',')

l.500 = c(as.vector(x.500$hgnc_symbol.x),as.vector(x.500$hgnc_symbol.y))
l.500 = paste(unique(l.500[which(l.500!='TP53')]),collapse = ',')

res = data.frame(SIGNATURE = c('STRING.ALL.TP53.X',"STRING.800.TP53.X",'STRING.TP53.500.X', 'STRING.ALL.TP53.AUTO'),GENES = c(l.all,l.800,l.500,l.all.auto.l))
res.tp53.auto = data.frame(SIGNATURE = c('STRING.ALL.TP53.AUTO'),GENES = c(l.all.auto.l))

res.auto = data.frame(GENES=l.all.auto)

write.table(res.auto,"~/Documents/PhD/GenderAnalysis/TP53_interactions/TP53.AUTO.Interactors.txt", sep = "\t", quote = F, row.names = F)

write.table(res,"~/Documents/PhD/GenderAnalysis/TP53_interactions/TP53.X.Signatures.txt",sep = ' ', quote = F, row.names = F)

write.table(res.tp53.auto ,"~/Documents/PhD/GenderAnalysis/TP53_interactions/TP53.AUTO.Signatures.txt",sep = ' ', quote = F, row.names = F)

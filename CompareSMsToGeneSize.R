library(dplyr)

output.dir = "~/Documents/PhD/GenderAnalysis/TCGA/Analysis/TCGAExpressionExplorerOutput/AllDataSets/"
sm.table = read.csv(paste(output.dir,'SM.Chi.Test.Propensity.csv',sep = ''))
sm.table = sm.table[which(sm.table$silent.mutations.females+sm.table$silent.mutations.males>=20),]
sm.table$total.silent.mut = sm.table$silent.mutations.females+sm.table$silent.mutations.males

gene.sizes = read.csv("~/Documents/PhD/GenderAnalysis/TCGA/Analysis/gene.sizes.csv")

gene.sizes = gene.sizes[,1:2]

colnames(gene.sizes) = c("Size",'feature')

gene.sizes %>% group_by(feature) %>% dplyr::summarise(SIZE = max(Size)) -> gene.sizes

sm.table = left_join(sm.table,gene.sizes)

sm.all = filter(sm.table,DATASET == 'ALL')

plot(sm.all$SIZE, sm.all$silent.mutations.females, main="Scatter", 
     xlab="Size ", ylab="SM ", pch=19)
abline(lm(sm.all$silent.mutations.females~sm.all$SIZE), col="red") # regression line (y~x) 
lines(lowess(sm.all$silent.mutations.females,sm.all$SIZE), col="blue") # lowess line (x,y)

cor(sm.all$SIZE, sm.all$silent.mutations.females,use = 'complete')
      
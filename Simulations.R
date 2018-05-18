rm(list = ls())

source('ExpressionPlotsFunctions.R')

###READ GENE lengths
gene.lengths = read.csv("~/Documents/PhD/GenderAnalysis/TCGA/Analysis/All.genes.lengths.csv")

###READ Dummy genes
dummy.genes = read.csv('~/Documents/PhD/GenderAnalysis/TCGA/Analysis/TCGAExpressionExplorerOutput/BLCA/DifferentialExpression/Mt.TP53.female-Mt.TP53.male.csv')

test.genes = dummy.genes[1:28,'X']

for (i in 1:1000){
  print(i)
  epa =GenerateRandomGeneSetBySize(test.genes,dummy.genes$X,data.base = gene.lengths)
  
}


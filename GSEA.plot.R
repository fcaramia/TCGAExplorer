rm (list = ls())
library(data.table)
library(ggplot2)
library(grid)

output.dir = "/home/fcaramia/Documents/PhD/GenderAnalysis/TCGA/Analysis/TCGAExpressionExplorerOutput/"
data = fread(paste(output.dir,'GSEA/gsea.pval.table.NEG.FC.csv',sep = ""))

data2 = data[which(data$contrast=='Normal.female-Wt.TP53.female'|data$contrast=='Normal.male-Wt.TP53.male')]

data2$gender = ifelse(data2$contrast=='Normal.female-Wt.TP53.female','Female','Male')
data2$stat = ifelse(data2$padj<=0.25,'Sig','NS')

grob <- grobTree(textGrob("Scatter plot", x=-1.9, hjust=0, y = 1,
                          gp=gpar(col="Blue", fontsize=13, fontface="italic")))


tiff(filename = paste(output.dir,'GSEA/gsea.plot.tiff',sep = ""),width = 600,height = 800, res = 100)
ggplot(data2, aes(x = NES, y = cancer, color = gender)) + geom_point(aes(size = -log10(pval),shape=stat)) + 
  xlim(-2,1) + ggtitle('Gene Set Enrichment, Negative Regulators p53') + scale_color_manual(values=alpha(c("hotpink1", "deepskyblue1"),1)) + 
  annotate("text", x = c(-1.7,0.8), y = c(11.4,11.4) , label = c("Wt-p53 Tumor",'Normal'), color = 'blue')
dev.off()

cosmic <- read.csv("~/Documents/PhD/Data/COSMIC/Census_allTue Nov 20 06_44_40 2018.csv")
cosmic <- cosmic[which(cosmic$Somatic=='yes'),]
write.csv(cosmic,"~/Documents/PhD/Data/COSMIC/Somatic_genes.csv", row.names = F)

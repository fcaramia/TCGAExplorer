cosmic <- read.csv("~/Documents/PhD/Data/COSMIC/Census_allTue Jan 17 04-15-44 2017.csv")
cosmic <- cosmic[which(substr(cosmic$Role.in.Cancer,1,4)=='onco'&cosmic$Somatic=='yes'),]
write.csv(cosmic,"~/Documents/PhD/Data/COSMIC/Somatic_oncogenes.csv", row.names = F)

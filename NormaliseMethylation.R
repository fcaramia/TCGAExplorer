rm(list = ls())
library(data.table)
library(dplyr)
library(FDb.InfiniumMethylation.hg19)
source("~/Documents/Work/analysis/RNA.Analysis/NormalisationFunctions.R")
source("~/Documents/Work/Projects/Sherene/Combined_Plotting/survival_functions.R")

output.dir = "~/Documents/PhD/GenderAnalysis/TCGA/Analysis/TCGAExpressionExplorerOutput/AllDataSets/Methylation/"
dir.create(output.dir,showWarnings = F)
escapers = c("PLCXD1","CSF2RA","IL3RA","SLC25A6","P2RY8","AKAP17A","TCONS_00017125","DHRSX","TCONS_00017281","CD99","XG","PRKX","LOC389906","HDHD1","TBL1X","MSL3","FRMPD4","TMSB4X","TRAPPC2","FANCB","CA5BP1","ZRSR2","AP1S2","SYAP1","TXLNG","PHEX","LOC100873065","EIF2S3","ZFX","IL1RAPL1","CXorf21","DMD","OTC","MED14","DDX3X","KDM6A","RBM3","USP27X","KDM5C","SMC1A","MAGED2","FAM104B","MTRNR2L10","LOC550643","FAAH2","EDA","RAB41","KIF4A","TEX11","FLJ44635","MAP2K4P1","XIST","JPX","FTX","TCONS_l2_00030263","BRWD3","KLHL4","DIAPH2","TSPAN6","TCONS_00017001","RBM41","COL4A6","ZCCHC16","AKAP14","LAMP2","C1GALT1C1","TCONS_00017461","XIAP","TENM1","TCONS_l2_00030350","MST4","SLC9A6","MAP7D3","TCONS_00017017","GABRA3")

methyl.x = fread('~/Documents/PhD/Data/TCGA_Xena/Methylation450/PANCAN_HumanMethyl_450.Xchrom') 
rownames(methyl.x) = methyl.x$sample
probenames = methyl.x$sample
methyl.x$sample = NULL
methyl.x = as.data.frame(methyl.x)

hm450 <- get450k()

probes <- hm450[probenames]
TSS = getNearestTSS(probes)
TSS$EscXInn = ifelse(TSS$nearestGeneSymbol%in%escapers, 'Escaper', 'Non-Escaper')
TSS$Esc = ifelse(TSS$nearestGeneSymbol%in%escapers, T, F)
gene.lenghts = read.csv('~/Documents/PhD/GenderAnalysis/TCGA/Analysis/gene.sizes.csv')

gene.lenghts %>% group_by(HGNC.symbol) %>% summarise(GeneLenght = max(Transcript.length..including.UTRs.and.CDS.)) -> gl 

TSS <- left_join(TSS, gl, by = c('nearestGeneSymbol'='HGNC.symbol'))

TSS$ProbePosPrctgSize = TSS$distance*100/TSS$GeneLenght

###


demo = fread("~/Documents/PhD/GenderAnalysis/TCGA/Analysis/all.TCGA.curated.clinical.csv")
demo$PATIENT_ID = toupper(demo$PATIENT_ID)

samples.ids = data.frame('SAMPLE_ID' = colnames(methyl.x))
samples.ids$PATIENT_ID = gsub("TCGA-([[:alnum:]]{2})-([[:alnum:]]{4}).*","\\2", samples.ids$SAMPLE_ID)
samples.ids = filter(samples.ids, PATIENT_ID%in%demo$PATIENT_ID)
samples.ids = left_join(samples.ids,demo[,c('CANCER_TYPE','PATIENT_ID','GENDER')])

samples.ids$SAMPLE_CODE = gsub("TCGA-([[:alnum:]]{2})-([[:alnum:]]{4}-([[:alnum:]]{2})).*","\\3", samples.ids$SAMPLE_ID)
samples.ids$SAMPLE_TYPE <- ifelse(as.integer(samples.ids$SAMPLE_CODE)>=20, 'Control',
                            ifelse(as.integer(samples.ids$SAMPLE_CODE)>=10, 'Normal', 'Tumor'))

###Use disparity set samples
datasets = c("ESCA","HNSC","LUSC","BLCA","LIHC","STAD","LGG","COAD","PAAD","READ","SKCM","LUAD")
demo = demo[demo$CANCER_TYPE%in%datasets,]
samples.ids = samples.ids[samples.ids$CANCER_TYPE%in%datasets,]
methyl.x = methyl.x[,as.vector(samples.ids$SAMPLE_ID)]
methyl.x = data.matrix(methyl.x)
##Check numbers
group_by(samples.ids,CANCER_TYPE,GENDER,SAMPLE_TYPE) %>% dplyr::summarise(SAMPLE_NUMBERS = n()) -> aux


##fill up the NAs

for (c in 1:dim(methyl.x)[1])
{
  methyl.x[c,which(is.na(methyl.x[c,]))] = median(methyl.x[c,],na.rm = T)
  
}

for (c in datasets[1]){
  #pdf(paste(output.dir,c,' RLE.pdf',sep = ""))
  samples.c = as.vector(samples.ids[samples.ids$CANCER_TYPE==c,'SAMPLE_ID'])
  samples.c.all = as.vector(samples.ids[samples.ids$CANCER_TYPE==c,])
  methyl.c = methyl.x[,samples.c]
  
  samples.c.male = samples.c.all[which(samples.c.all$GENDER=='male'),]
  samples.c.female = samples.c.all[which(samples.c.all$GENDER=='female'),]
 #DoMultPCAPlots(data = t(methyl.c),annot = samples.c.all,dir = output.dir,title = paste(c,'PCA'),color.cols = c('GENDER'),all.samps.by = c('GENDER'), filename = c)
  DoRLEPlot(data = methyl.c[,as.vector(samples.c.female$SAMPLE_ID)], annot = samples.c.female,color.col = 'SAMPLE_TYPE', title = paste(c,'females','RLE') )
  DoRLEPlot(data = methyl.c[,as.vector(samples.c.male$SAMPLE_ID)], annot = samples.c.male,color.col = 'SAMPLE_TYPE', title = paste(c,'males','RLE') )
  
  DoRLEPlot(data = 2*abs(methyl.c[,as.vector(samples.c.female$SAMPLE_ID)]-0.5), annot = samples.c.female,color.col = 'SAMPLE_TYPE', title = paste(c,'females transformed ','RLE') )
  
  
  
  fem.methyl = methyl.c[,as.vector(samples.c.female$SAMPLE_ID)]
  
  probesMeds = rowMedians(fem.methyl, na.rm = T)
  
  TSS$probesMeds = probesMeds
  
  TSS.s = filter(TSS, ProbePosPrctgSize  < 100)
  
  plot(TSS.s$ProbePosPrctgSize, TSS.s$probesMeds, main="Distance to TSS by % vs Methylation", 
       xlab="Distance to TSS ", ylab="B value ")
  
  cor(TSS.s$distance, TSS.s$probesMeds)
  
  ggplot(data = TSS,aes(x = 2*abs(TSS$probesMeds-0.5) , colour=TSS$EscXInn)) + geom_density() +
    scale_color_brewer(palette = 'Set1') + ggtitle(label = 'Escapers v non escapers') + theme_minimal()
  
  ggplot(data = TSS,aes(x =TSS$probesMeds , colour=TSS$EscXInn)) + geom_density() +
    scale_color_brewer(palette = 'Set1') + ggtitle(label = 'Escapers v non escapers') + theme_minimal()
  
  
  TSS$prox = ifelse(TSS$ProbePosPrctgSize<5, 'TSS', 'Body')
  TSS$disp = paste(TSS$EscXInn, TSS$prox)
  
  TSS.d = TSS[!is.na(TSS$ProbePosPrctgSize),]
  TSS.d$probesMeds.transformed = 2*abs(TSS.d$probesMeds - 0.5)
  
  bplot(mat = TSS.d,datcol = 'probesMeds', typecol = 'disp', title = 'Meth')
  
  #dev.off()
}






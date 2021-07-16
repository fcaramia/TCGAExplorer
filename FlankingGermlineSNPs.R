rm(list = ls())
library(data.table)
library(dplyr)
source("ExpressionPlotsFunctions.R")

escapers = c("PLCXD1","CSF2RA","IL3RA","SLC25A6","P2RY8","AKAP17A","TCONS_00017125","DHRSX","TCONS_00017281","CD99","XG","PRKX","LOC389906","HDHD1","TBL1X","MSL3","FRMPD4","TMSB4X","TRAPPC2","FANCB","CA5BP1","ZRSR2","AP1S2","SYAP1","TXLNG","PHEX","LOC100873065","EIF2S3","ZFX","IL1RAPL1","CXorf21","DMD","OTC","MED14","DDX3X","KDM6A","RBM3","USP27X","KDM5C","SMC1A","MAGED2","FAM104B","MTRNR2L10","LOC550643","FAAH2","EDA","RAB41","KIF4A","TEX11","FLJ44635","MAP2K4P1","XIST","JPX","FTX","TCONS_l2_00030263","BRWD3","KLHL4","DIAPH2","TSPAN6","TCONS_00017001","RBM41","COL4A6","ZCCHC16","AKAP14","LAMP2","C1GALT1C1","TCONS_00017461","XIAP","TENM1","TCONS_l2_00030350","MST4","SLC9A6","MAP7D3","TCONS_00017017","GABRA3")

output.dir = "~/Documents/PhD/GenderAnalysis/TCGA/Analysis/TCGAExpressionExplorerOutput/AllDataSets/GermlineRMAF/"
dir.create(output.dir,showWarnings = F)

purity = fread("~/Documents/PhD/Data/TCGA_2016_01_28_BROAD/Purity.Data/Purity_Ploidy_All_Samples_4-17-15.csv")

tumor.snps = fread("~/Documents/PhD/GenderAnalysis/TCGA/Analysis/Germline.Tumor.ExMut.v2.csv")
tumor.snps = tumor.snps[tumor.snps$CHROMOSOME=='X'&nchar(tumor.snps$REFERENCE_ALLELE)<6&nchar(tumor.snps$TUMOR_SEQ_ALLELE2)<6,]

normals.snps = fread("~/Documents/PhD/GenderAnalysis/TCGA/Analysis/Germline.Normals.ExMut.V2.csv")
normals.snps = normals.snps[normals.snps$CHROMOSOME=='X'&nchar(normals.snps$REFERENCE_ALLELE)<6&nchar(normals.snps$TUMOR_SEQ_ALLELE2)<6,]

clinical.data = read.csv("~/Documents/PhD/GenderAnalysis/TCGA/Analysis/all.TCGA.curated.clinical.csv", as.is = T)
clinical.data$PATIENT_ID = toupper(clinical.data$PATIENT_ID)

normals.snps$SNP.ID = paste(normals.snps$PATIENT_ID,normals.snps$CHROMOSOME,normals.snps$START_POSITION,normals.snps$REFERENCE_ALLELE, sep = '.')
tumor.snps$SNP.ID = paste(tumor.snps$PATIENT_ID,tumor.snps$CHROMOSOME,tumor.snps$START_POSITION,tumor.snps$REFERENCE_ALLELE, sep = '.')

#####################################
###Duplicates##########
normals.snps$NORMAL_SUM = normals.snps$REF_COUNT + normals.snps$ALT_COUNT
normals.snps = normals.snps[!is.na(normals.snps$NORMAL_SUM),]
tumor.snps$TUMOR_SUM = tumor.snps$REF_COUNT + tumor.snps$ALT_COUNT
tumor.snps = tumor.snps[!is.na(tumor.snps$TUMOR_SUM),]

normals.snps %>% group_by(SNP.ID) %>% dplyr::summarise(s=max(NORMAL_SUM)) %>% 
  dplyr::mutate(U_ID = paste(SNP.ID,s,sep = '.')) -> aux

normals.snps$U_ID = paste(normals.snps$SNP.ID, normals.snps$NORMAL_SUM,sep = '.')
normals.snps %>% filter(U_ID %in% aux$U_ID) -> normals.snps

tumor.snps %>% group_by(SNP.ID) %>% dplyr::summarise(s=max(TUMOR_SUM)) %>% 
  dplyr::mutate(U_ID = paste(SNP.ID,s,sep = '.')) -> aux

tumor.snps$U_ID = paste(tumor.snps$SNP.ID, tumor.snps$TUMOR_SUM,sep = '.')
tumor.snps %>% filter(U_ID %in% aux$U_ID) -> tumor.snps

normals.snps = normals.snps[!duplicated(normals.snps$U_ID),]
tumor.snps = tumor.snps[!duplicated(tumor.snps$U_ID),]
################################


normals.snps %>% mutate(RMAF_NORMALS = ifelse(NORMAL_SUM>=10,ALT_COUNT/NORMAL_SUM,0.0)) -> normals.snps
tumor.snps %>% mutate(RMAF_TUMORS = ifelse(TUMOR_SUM>=10,ALT_COUNT/TUMOR_SUM,0.0)) -> tumor.snps

tumor.snps %>% left_join(normals.snps[,c('SNP.ID','NORMAL_SUM','RMAF_NORMALS')]) -> snps

snps$RMAF_DELTA = snps$RMAF_NORMALS - snps$RMAF_TUMORS

snps%>% left_join(clinical.data[,c('PATIENT_ID','GENDER','CANCER_TYPE')]) -> snps

snps$ABS_DELTA = abs(snps$RMAF_DELTA)

snps = snps[snps$NORMAL_SUM>=10&snps$TUMOR_SUM>=10,]

snps$individual_id = gsub("TCGA\\-([[:alnum:]]{2})\\-([[:alnum:]]{4})\\-([[:alnum:]]{2}).*","TCGA\\-\\1\\-\\2", snps$SAMPLE_ID)

snps%>% left_join(purity) -> snps

snps$PUR = ifelse(snps$purity>.75,'high',ifelse(snps$purity<.5,'low','med'))

####Violin plot
pdf(paste(output.dir,'ViolinMvFHom.pdf',sep = ''))
DoMutNumViolinPlotGender(data = snps,data.col = 'RMAF_DELTA',fill.col = 'GENDER', 
                         title.txt = 'Delta', x.lab = 'Gender', y.lab = 'delta.abs',show.points = T,points.alpha = .05)
dev.off()

##Scatter plot
pdf(paste(output.dir,'ScatterMvFHom.pdf',sep = ''))
b <- ggplot(snps, aes(x = RMAF_TUMORS, y = RMAF_NORMALS))
b + geom_point(aes(color = GENDER))
dev.off()
###Genome Wide Plot


library(biovizBase)
XchromCyto <- getIdeogram("hg19", cytobands = TRUE,subchr = 'chrX')
Xdf = data.frame(XchromCyto, stringsAsFactors = F)
n = length(levels(as.factor(as.vector(Xdf$name))))
cols = rainbow(n, s=.6, v=.9)[sample(1:n,n)]
Xdf$color = cols
Xdf$offset = rep(c(0,-0.22),20)

pdf(paste(output.dir,'XChromMvFHom.pdf',sep = ''),width = 12,height = 9)
p <- ggplot(data = snps,aes(x=START_POSITION,y = ABS_DELTA)) 
p +  geom_point(size = 1.5, alpha=0.05) + facet_grid(CANCER_TYPE~GENDER,scales =  'free_x') + 
  geom_segment(data = Xdf,inherit.aes = F,mapping = aes(x=start,y=-.40, xend = end, yend=-.40,color=name), size=4, show.legend = F) +
  scale_color_manual(values = Xdf$color) +scale_y_continuous(limits = c(-.70,1.15)) + 
  geom_text(data = Xdf,inherit.aes = F,show.legend = F,mapping = aes(x=start,y = -.4+offset,label=name),size = 1.75,hjust=-0.1)

dev.off()

# Call the palette with a number
pdf(paste(output.dir,'Density2DMvFHom.pdf',sep = ''),width = 12,height = 12)
ggplot(snps, aes(x=RMAF_TUMORS, y=RMAF_NORMALS) ) +
  stat_density_2d(aes(fill = ..level..), geom = "polygon") + scale_fill_viridis_c() + facet_grid(.~GENDER)
dev.off()

#annotate genes

# Load the library

snps.fem = snps[snps$GENDER=='female',]
snps.fem = snps.fem[order(snps.fem$START_POSITION),]

library(biomaRt)

mart <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl",host = 'asia.ensembl.org')
genes<-getBM(c("hgnc_symbol","ensembl_gene_id","chromosome_name","start_position","end_position"),filters = c('chromosome_name'),values = c('X') , mart=mart)
genes = genes[order(genes$start_position),]
library(GenomicRanges)

ref <- GRanges(seqnames = Rle(genes$chromosome_name), ranges = IRanges(start = genes$start_position, end = genes$end_position), names = genes$hgnc_symbol)
snps.gr <- GRanges(seqnames = Rle(snps.fem$CHROMOSOME),ranges = IRanges(start = snps.fem$START_POSITION, end = snps.fem$START_POSITION+1))
olap <- findOverlaps(snps.gr, ref,type = 'within',select = 'first')

snps.fem = snps.fem[!is.na(olap),]
olap = olap[!is.na(olap)]
snps.fem$Gene = ref[olap]$names

snps.fem$ESCAPE = ifelse(snps.fem$Gene%in%escapers,'Yes','No')

####Do Plots

####Violin plot
pdf(paste(output.dir,'ViolinESCAPEvNoESCAPE.pdf',sep = ''))
DoMutNumViolinPlotGender(data = snps.fem,data.col = 'RMAF_DELTA',fill.col = 'ESCAPE', 
                         title.txt = 'Delta', x.lab = 'Gender', y.lab = 'delta.abs',show.points = T,points.alpha = .05)
dev.off()

##Scatter plot
pdf(paste(output.dir,'ScatterESCAPEvNoESCAPE.pdf',sep = ''))
b <- ggplot(snps.fem, aes(x = RMAF_TUMORS, y = RMAF_NORMALS))
b + geom_point(aes(color = ESCAPE))
dev.off()
###Genome Wide Plot


library(biovizBase)
XchromCyto <- getIdeogram("hg19", cytobands = TRUE,subchr = 'chrX')
Xdf = data.frame(XchromCyto, stringsAsFactors = F)
n = length(levels(as.factor(as.vector(Xdf$name))))
cols = rainbow(n, s=.6, v=.9)[sample(1:n,n)]
Xdf$color = cols
Xdf$offset = rep(c(0,-0.22),20)

pdf(paste(output.dir,'XChromESCAPEvNoESCAPE.pdf',sep = ''),width = 12,height = 9)
p <- ggplot(data = snps.fem,aes(x=START_POSITION,y = ABS_DELTA)) 
p +  geom_point(size = 1.5, alpha=0.05) + facet_grid(CANCER_TYPE~ESCAPE,scales =  'free_x') + 
  geom_segment(data = Xdf,inherit.aes = F,mapping = aes(x=start,y=-.40, xend = end, yend=-.40,color=name), size=4, show.legend = F) +
  scale_color_manual(values = Xdf$color) +scale_y_continuous(limits = c(-.70,1.15)) + 
  geom_text(data = Xdf,inherit.aes = F,show.legend = F,mapping = aes(x=start,y = -.4+offset,label=name),size = 1.75,hjust=-0.1)

dev.off()

# Call the palette with a number
pdf(paste(output.dir,'Density2DESCAPEvNoESCAPE.pdf',sep = ''),width = 12,height = 12)
ggplot(snps.fem, aes(x=RMAF_TUMORS, y=RMAF_NORMALS) ) +
  stat_density_2d(aes(fill = ..level..), geom = "polygon") + scale_fill_viridis_c() + facet_grid(.~ESCAPE)
dev.off()

pdf(paste(output.dir,'Density2DPurity.pdf',sep = ''),width = 12,height = 12)
ggplot(snps.fem[which(!is.na(snps.fem$PUR)),], aes(x=RMAF_TUMORS, y=RMAF_NORMALS) ) +
  stat_density_2d(aes(fill = ..level..), geom = "polygon") + scale_fill_viridis_c() + facet_grid(.~PUR)
dev.off()

snps.fem$PUR.ESC = paste(snps.fem$PUR,snps.fem$ESCAPE,sep = '.')
pdf(paste(output.dir,'Density2DPurityESC.pdf',sep = ''),width = 12,height = 12)
ggplot(snps.fem[which(!is.na(snps.fem$PUR)),], aes(x=RMAF_TUMORS, y=RMAF_NORMALS) ) +
  stat_density_2d(aes(fill = ..level..), geom = "polygon") + scale_fill_viridis_c() + facet_grid(.~PUR.ESC)
dev.off()



#####p53 interactors
p53.int = read.csv("~/Documents/PhD/GenderAnalysis/TP53_interactions/Candidates.2018.03.01.csv")

p53.int = p53.int$ID

snps.fem$TP53.INT = ifelse(snps.fem$Gene%in%p53.int,'Yes','No')
#####PLOTS
####Violin plot
pdf(paste(output.dir,'ViolinP53Int.pdf',sep = ''))
DoMutNumViolinPlotGender(data = snps.fem,data.col = 'RMAF_DELTA',fill.col = 'TP53.INT', 
                         title.txt = 'Delta', x.lab = 'Gender', y.lab = 'delta.abs',show.points = T,points.alpha = .05)
dev.off()

##Scatter plot
pdf(paste(output.dir,'ScatterP53Int.pdf',sep = ''))
b <- ggplot(snps.fem, aes(x = RMAF_TUMORS, y = RMAF_NORMALS))
b + geom_point(aes(color = TP53.INT))
dev.off()
###Genome Wide Plot

# Call the palette with a number
pdf(paste(output.dir,'Density2DP53.INT.pdf',sep = ''),width = 12,height = 12)
ggplot(snps.fem, aes(x=RMAF_TUMORS, y=RMAF_NORMALS) ) +
  stat_density_2d(aes(fill = ..level..), geom = "polygon") + scale_fill_viridis_c() + facet_grid(.~TP53.INT)
dev.off()

####Mutant p53 samples

###Mutations####
mutations = fread("~/Documents/PhD/GenderAnalysis/TCGA/Analysis/full.reduced.all.raw.TCGA.curated.mutations.csv")
mutations$PATIENT_ID = gsub("TCGA\\-([[:alnum:]]{2})\\-([[:alnum:]]{4}).*","\\2", mutations$SAMPLE_ID)
###############

filter.out = c("Silent",'Intron','IGR',"In_Frame_Ins" ,"In_Frame_Del", "lincRNA" )
mutations %>% 
  filter(!(VARIANT_CLASSIFICATION%in%filter.out),HUGO_SYMBOL=='TP53') -> mutations

snps.fem$TP53.MUT = ifelse(snps.fem$PATIENT_ID%in%mutations$PATIENT_ID,'Mt','Wt')

#####PLOTS
####Violin plot
pdf(paste(output.dir,'ViolinP53Status.pdf',sep = ''))
DoMutNumViolinPlotGender(data = snps.fem,data.col = 'RMAF_DELTA',fill.col = 'TP53.MUT', 
                         title.txt = 'Delta', x.lab = 'Gender', y.lab = 'delta.abs',show.points = T,points.alpha = .05)
dev.off()

##Scatter plot
pdf(paste(output.dir,'ScatterP53Status.pdf',sep = ''))
b <- ggplot(snps.fem, aes(x = RMAF_TUMORS, y = RMAF_NORMALS))
b + geom_point(aes(color = TP53.MUT))
dev.off()
###Genome Wide Plot

# Call the palette with a number
pdf(paste(output.dir,'Density2DP53.STATUS.pdf',sep = ''),width = 12,height = 12)
ggplot(snps.fem, aes(x=RMAF_TUMORS, y=RMAF_NORMALS) ) +
  stat_density_2d(aes(fill = ..level..), geom = "polygon") + scale_fill_viridis_c() + facet_grid(.~TP53.MUT)
dev.off()



rm (list = ls())
library(data.table)
library(reshape)
library(ggplot2)
library(dplyr)
library(survival)
source("~/Documents/Work/analysis/RNA.Analysis/NormalisationFunctions.R")
###Read Achilles data####
ataris = fread("~/Documents/PhD/Data/CCLE/Achilles/Achilles_QC_v2.4.3.rnai.Gs.gct")
dat = melt(ataris)
dat$TISSUE = as.factor(gsub(pattern = ".*_(.*)",x = dat$variable,replacement = '\\1'))
dat = dat[dat$TISSUE%in%c('BREAST','OVARY','KIDNEY','LIVER',"LUNG","OESOPHAGUS",'PANCREAS','SKIN',"STOMACH", 'INTESTINE','PLEURA', 'TRACT')]



####READ CCLE DATA####
ccle.muts = fread("~/Documents/PhD/Data/CCLE/CCLE_DepMap_18q3_maf_20180718.txt")
ccle.muts$GENE.CL = paste(ccle.muts$Hugo_Symbol,ccle.muts$Tumor_Sample_Barcode,sep = '.')
ccle.annot = fread("~/Documents/PhD/Data/CCLE/CCLE_sample_info_file_2012-10-18.txt")
P.Coding.X.genes = read.csv("~/Documents/PhD/GenderAnalysis/Genes.Coding.X.csv")
p53.muts = ccle.muts[ccle.muts$Hugo_Symbol=='TP53',]


####P53 interactors
p53.int = read.csv("~/Documents/PhD/GenderAnalysis/TP53_interactions/Candidates.2018.12.20.full.csv")
#p53.high = as.vector(p53.int[which(p53.int$score>=0.6),'ID'])
#p53.med = as.vector(p53.int[which(p53.int$score>=0.4),'ID'])
p53.low = as.vector(p53.int[which(p53.int$experimentally_determined_interaction>=0.3),'GeneSymbol'])
p53.low = as.vector(p53.int$GeneSymbol)

###Read RNA SEQ
ccle.rna.norm = as.data.frame(fread("~/Documents/PhD/Data/CCLE/CCLE.NORM.COUNTS.csv"))
rownames(ccle.rna.norm) = ccle.rna.norm$Gene.ID
ccle.rna.norm$Gene.ID = NULL
###Keep useful genes and cell lines###
keep.cels = which(colnames(ccle.rna.norm)%in%dat$variable)
keep.genes = which(rownames(ccle.rna.norm)%in%dat$Description)

ccle.rna = ccle.rna.norm[keep.genes,keep.cels]

####Gene annotation for EntrezGeneID
gene.annot = read.csv("~/Documents/PhD/GenderAnalysis/TCGA/Analysis/annotation.all.genes.csv",as.is = T)
gene.annot = gene.annot[,-1]
gene.annot = gene.annot[gene.annot$Chrm!='',]
gene.annot = gene.annot[order(gene.annot$Chrm,gene.annot$Entrez.Gene.ID),]

p53.int$EntrezID = gene.annot[match(p53.int$GeneSymbol,gene.annot$Approved.Symbol),'Entrez.Gene.ID']


###Compute signature 
sig = p53.int[which(p53.int$GeneSymbol%in%rownames(ccle.rna)&p53.int$GeneSymbol%in%p53.low),]
colnames(sig) = c('probe','score','EntrezGene.ID')
sig$coefficient = 1

score = sig.score(x = sig,data = data.matrix(t(ccle.rna)),annot = sig)
score = as.data.frame(score$score)
score$CL = rownames(score)
score$score = score$`score$score`

ccle.rna$Gene.ID = rownames(ccle.rna)
dat.ccle = melt(ccle.rna)


###Merge achilles and p53 rna sig
dat$P53.STRING.SIG = score[match(dat$variable,score$CL),'score']

###Merge achilles and ccle rnaseq 
dat.ccle$Gene.CL = paste(dat.ccle$Gene.ID,dat.ccle$variable,sep = '_')
dat$Gene.CL = paste(dat$Description,dat$variable,sep = '_')
dat$CCLE = dat.ccle[match(dat$Gene.CL,dat.ccle$Gene.CL),'value']


####Check siRNA, RNASeq correlation

# ggplot(dat[which(abs(dat$value)>2.5),], aes(x=value, y=CCLE)) +
#   geom_point(size=0.5, shape=23)

# cor.test(~ CCLE+abs(value), data = dat[which(abs(dat$value)>2.5),])

###The majority of genes with high abs(siRNA) have high RNA Seq 
###Use cutoffs of abs(siRNA) > 2.5 and RNA Seq > 0
###Use cutoff of RNASeq > 0 if not enough genes

###Merge with cosmic#####
cosmic <- read.csv("~/Documents/PhD/Data/COSMIC/Census_allTue Nov 20 06_44_40 2018.csv")
dat$ROLE = cosmic[match(dat$Description,cosmic$Gene.Symbol),'Role.in.Cancer']
dat[which(dat$ROLE==''),'ROLE'] = NA
dat$TISSUE.COSMIC = cosmic[match(dat$Description,cosmic$Gene.Symbol),'Tumour.Types.Somatic.']
dat$GENE.CL = paste(dat$Description,dat$variable,sep='.')
dat$GENE.STATUS = as.factor(ifelse(dat$GENE.CL%in%ccle.muts$GENE.CL,'Mt','Wt'))
dat$COSMIC = as.factor(ifelse(dat$Description%in%cosmic$Gene.Symbol,'Cosmic','No'))

dat$TISSUE.STATUS = as.factor(paste(dat$TISSUE,dat$GENE.STATUS,sep='.'))

dat$VAR.CLASS = ccle.muts[match(dat$GENE.CL,ccle.muts$GENE.CL),'Variant_Classification']
dat[which(is.na(dat$VAR.CLASS)),'VAR.CLASS'] = 'Wt'

dat$IS.DEL = ccle.muts[match(dat$GENE.CL,ccle.muts$GENE.CL),'isDeleterious']
dat$IS.DEL = as.factor(ifelse(dat$IS.DEL==T,'Mt-DEL',ifelse(is.na(dat$IS.DEL),'Wt','Mt')))


####Check COSMIC reported in tissue

dat$COSMIC.IN.TISSUE =  ifelse(dat$TISSUE=='KIDNEY',grepl('kidney', dat$TISSUE.COSMIC)|grepl('renal', dat$TISSUE.COSMIC),
                        ifelse(dat$TISSUE=='LIVER',grepl('liver', dat$TISSUE.COSMIC)|grepl('hepatocellular', dat$TISSUE.COSMIC),
                        ifelse(dat$TISSUE=='LUNG',grepl('lung', dat$TISSUE.COSMIC),
                        ifelse(dat$TISSUE=='BREAST',grepl('breast', dat$TISSUE.COSMIC),
                        ifelse(dat$TISSUE=='OVARY',grepl('ovar', dat$TISSUE.COSMIC),
                        ifelse(dat$TISSUE=='OESOPHAGUS',grepl('oesophageal', dat$TISSUE.COSMIC),
                        ifelse(dat$TISSUE=='PANCREAS',grepl('pancrea', dat$TISSUE.COSMIC),
                        ifelse(dat$TISSUE=='SKIN',grepl('skin', dat$TISSUE.COSMIC)|grepl('melanoma', dat$TISSUE.COSMIC),
                        ifelse(dat$TISSUE=='STOMACH',grepl('stomach', dat$TISSUE.COSMIC),
                        ifelse(dat$TISSUE=='INTESTINE',grepl('colorectal', dat$TISSUE.COSMIC)|grepl('intestine', dat$TISSUE.COSMIC)|grepl('intestinal', dat$TISSUE.COSMIC),
                        ifelse(dat$TISSUE=='PLEURA',grepl('lung', dat$TISSUE.COSMIC),
                        ifelse(dat$TISSUE=='TRACT',grepl('bladder', dat$TISSUE.COSMIC)|grepl('biliary', dat$TISSUE.COSMIC),
                        ifelse(grepl('other', dat$TISSUE.COSMIC),T,F)))))))))))))
                                                                                
                                                                  


# ###Explore siRNA based on mutations of Genes and role in cancer
# dat.plot = dat[which(abs(dat$value)>=2&dat$CCLE>=0),]
# p <- ggplot(data = dat.plot, aes(y = value , x = COSMIC.IN.TISSUE  , color = GENE.STATUS, alpha= 0.2))
# p + geom_boxplot()  + facet_wrap(dat.plot$ROLE)





####Read p53 String SET



#dat$P53.HIGH = ifelse(dat$Description%in%p53.high,'p53-int','No p53-int')
#dat$P53.MED = ifelse(dat$Description%in%p53.med,'p53-int','No p53-int')
dat$P53.LOW = as.factor(ifelse(dat$Description%in%p53.low,T,F))
dat$P53.MUT = as.factor(ifelse(dat$variable%in%p53.muts$Tumor_Sample_Barcode,'mt-p53','wt-p53'))
dat$X = ifelse(dat$Description%in%P.Coding.X.genes$hgnc_symbol,'X','AUTOSOME')

dat$COMB = as.factor(paste(dat$P53.LOW,dat$P53.MUT,sep = ' | '))
dat = dat[!is.na(dat$value),]

dat$INT.ROLE = as.factor(paste(dat$P53.LOW,dat$ROLE,sep='.'))
dat$EXPRSS = as.factor(ifelse(dat$CCLE>=median(dat$CCLE,na.rm = T),'High','Low'))
dat$MUT.EXPRESS = as.factor(paste(dat$GENE.STATUS,dat$EXPRSS,sep = '-'))
###### Plot Scores ####

dat$p53INT.MUT = dat$GENE.STATUS=='Mt'&dat$P53.LOW==T
dat$COSMIC.MUT = dat$COSMIC=='Cosmic'&dat$GENE.STATUS=='Mt'
dat$p53INT.MUT.CL = dat$GENE.STATUS=='Mt'&dat$P53.LOW==T&dat$COSMIC.IN.TISSUE
dat$COSMIC.MUT.CL = dat$COSMIC=='Cosmic'&dat$GENE.STATUS=='Mt'& dat$COSMIC.IN.TISSUE
dat$GENE.STATUS.P53 = paste(dat$GENE.STATUS,'-gene','.',dat$P53.MUT,sep = '')

###count mutations of p53 string genes

ccle.muts %>% filter(Hugo_Symbol%in%p53.low) -> p53.string.muts
p53.string.muts %>% group_by(Tumor_Sample_Barcode) %>% dplyr::summarise(N.P = n()) -> p53.string.muts.summ
dat$p53.string.em.muts = p53.string.muts.summ[match(dat$variable,p53.string.muts.summ$Tumor_Sample_Barcode),'N.P']
dat[which(is.na(dat$p53.string.em.muts)),'p53.string.em.muts'] = 0


q.sig = quantile(dat$P53.STRING.SIG,probs = c(.45,.55),na.rm = T)
#q.sig = quantile(dat$P53.STRING.SIG,probs = c(.50),na.rm = T)

dat$P53.STRING.SIG.d = ifelse(dat$P53.STRING.SIG<=q.sig[1],'Low',ifelse(dat$P53.STRING.SIG>=q.sig[2],'High','Med'))
#dat$P53.STRING.SIG.d = ifelse(dat$P53.STRING.SIG<q.sig[1],'Low','High')
dat[is.na(dat$P53.STRING.SIG),'P53.STRING.SIG.d'] = NA

###Explore siRNA based on mutations of Genes and role in cancer with p53 interactors
# dat.plot = dat[which(!is.na(dat$CCLE)&abs(dat$value)>=0&dat$CCLE>=0&dat$P53.LOW==T),]
# p <- ggplot(data = dat.plot, aes(y = value , x = GENE.STATUS.P53, color = GENE.STATUS.P53))
# p + geom_boxplot() + facet_wrap(~dat.plot$Description) + theme(axis.text.x=element_text(angle=90, hjust=1))
# 
# dat.plot = dat[which(!is.na(dat$CCLE)&abs(dat$value)>=0&dat$CCLE>=0&dat$P53.LOW==T),]
# p <- ggplot(data = dat.plot, aes(y = value , x = CCLE, color = GENE.STATUS.P53))
# p + geom_point()
# 
# dat.plot = dat[which(!is.na(dat$CCLE)&abs(dat$value)>=0&dat$CCLE>=0&dat$P53.LOW==T),]
# p <- ggplot(data = dat.plot, aes(y = value , x = CCLE, color = GENE.STATUS))
# p + geom_boxplot()


####Test muts in p53.string set and p53 KD
dat.plot = dat[which(!is.na(dat$CCLE)&abs(dat$value)>=0&dat$CCLE>=0&dat$Description=='TP53'),]
dat.plot$p53.string.em.muts.d = ifelse(dat.plot$p53.string.em.muts==0,'Low',
                                       ifelse(dat.plot$p53.string.em.muts<10,'Low','High'))
p <- ggplot(data = dat.plot, aes(y = value , fill = p53.string.em.muts.d))
p + geom_boxplot() + facet_wrap(~dat.plot$P53.MUT) + guides(fill=guide_legend(title="Mutation incidence\n in p53 String set")) +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) + 
        labs(y = 'Cancer Growth and Development' )

test.dat = dat.plot[which(dat.plot$P53.MUT=='wt-p53'),c('value','p53.string.em.muts.d')]
test.dat$p53.string.em.muts.d = as.factor(ifelse(test.dat$p53.string.em.muts.d=='High',1,0))
wilcox.test(value~p53.string.em.muts.d,data=data.matrix(test.dat))



dat.plot = dat[which(!is.na(dat$CCLE)&abs(dat$value)>=0&dat$CCLE>=0&dat$Description=='TP53'&dat$P53.STRING.SIG.d!='Med'),]
p <- ggplot(data = dat.plot, aes(y = value , fill = P53.STRING.SIG.d))
p + geom_boxplot() + facet_wrap(~dat.plot$P53.MUT) 


test.dat = dat.plot[which(dat.plot$P53.MUT=='wt-p53'),c('value','P53.STRING.SIG.d')]
test.dat$P53.STRING.SIG.d = as.factor(ifelse(test.dat$P53.STRING.SIG.d=='High',1,0))
wilcox.test(value~P53.STRING.SIG.d,data=data.matrix(test.dat))

as.factor(dat[which(dat$p53INT.MUT==T),'Description'])

dat.plot = dat[which(!is.na(dat$CCLE)&abs(dat$value)>=0&dat$CCLE>=1),]
p <- ggplot(data = dat.plot, aes(y = value , x = p53INT.MUT.CL,color=p53INT.MUT.CL))
p + geom_boxplot()



summary(dat[which(dat$P53.LOW==T&!is.na(dat$TISSUE.COSMIC)),'TISSUE.COSMIC'])


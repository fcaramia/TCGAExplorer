rm(list = ls())
library(dplyr)
library(ggrepel)
source('ExpressionPlotsFunctions.R')
source('PropensityScoresFunctions.R')
library(gtools)
output.dir = "~/Documents/PhD/GenderAnalysis/TCGA/Analysis/TCGAExpressionExplorerOutput/AllDataSets/"




#Datasets and directories and variables
datasets = c('ALL',"ESCA","HNSC","LUSC","BLCA","LIHC","STAD","LGG","COAD","PAAD","READ","SKCM")
#datasets = c("READ")

RMAF = read.csv(paste(output.dir,'Computed.RMAF.csv',sep = ''))

#Clean low reads
RMAF$SUM = RMAF$REF_COUNT + RMAF$ALT_COUNT

RMAF %>% mutate(RMAF = ifelse(SUM>=10,ALT_COUNT/SUM,0.0)) %>% filter(!(EXPRS.GENE=='NOT.EXPRS'&SUM<10)) %>%
  filter(!(SUM<10&Z.SCORE>-4))-> RMAF

RMAF[which(RMAF$Z.SCORE<=-4&RMAF$SUM<10),'RMAF'] = 1.0
RMAF$LOGIT.RMAF = logit(x = RMAF$RMAF,min = 0-0.1 , max = 1+0.1)
#RMAF[which(RMAF$Z.SCORE<=-4&RMAF$SUM<5),'RMAF'] = 1.0
#DoDensityPlotGenderRMAF(dat = RMAF,dat.col = 'RMAF',fill.col = 'GENDER',title.txt = 'RMAF')
#DoDensityPlotGenderRMAF(dat = RMAF,dat.col = 'Z.SCORE',fill.col = 'GENDER',title.txt = 'Z-Scores')
#Define Expressed Mutation RMAF
RMAF$EM.RMAF = ifelse(RMAF$RMAF>=0.75,1,0)
RMAF$SM.RMAF = ifelse(RMAF$RMAF<=0.20,1,0)


RMAF %>% group_by(GENDER,EM.RMAF) %>% dplyr::summarise(N.EM=n()) -> EM.NUMS
RMAF %>% group_by(GENDER,SM.RMAF) %>% dplyr::summarise(N.SM=n()) -> SM.NUMS

full_join(EM.NUMS,SM.NUMS) -> TOT.NUMS

write.csv(TOT.NUMS,paste(output.dir,'Total.SM.EM.by.gender.csv',sep = ''), row.names = F)


#Use mutations present in 10+ patients
RMAF %>% group_by(HUGO_SYMBOL, PATIENT_ID) %>% dplyr::summarise(s=max(RMAF)) %>% 
  dplyr::mutate(U_ID = paste(HUGO_SYMBOL, PATIENT_ID,s,sep = '.')) -> aux

RMAF$U_ID = paste(RMAF$HUGO_SYMBOL, RMAF$PATIENT_ID,RMAF$RMAF,sep = '.')
RMAF %>% filter(U_ID %in% aux$U_ID) -> RMAF

RMAF %>% group_by(HUGO_SYMBOL, PATIENT_ID) %>% dplyr::summarise(s=max(START_POSITION)) %>% 
  dplyr::mutate(U_ID = paste(HUGO_SYMBOL, PATIENT_ID,s,sep = '.')) -> aux

RMAF$U_ID = paste(RMAF$HUGO_SYMBOL, RMAF$PATIENT_ID,RMAF$START_POSITION,sep = '.')
RMAF %>% filter(U_ID %in% aux$U_ID) -> RMAF

RMAF %>% group_by(HUGO_SYMBOL) %>% dplyr::summarise(n=n()) %>% filter(n>=10) -> aux
RMAF %>% filter(HUGO_SYMBOL %in% aux$HUGO_SYMBOL) -> RMAF
rm(aux)

DNA.VAF = read.csv("~/Documents/PhD/GenderAnalysis/TCGA/Analysis/full.reduced.all.raw.TCGA.curated.mutations.csv", as.is=T)
DNA.VAF %>% filter(CANCER_TYPE%in%datasets) -> DNA.VAF
DNA.VAF %>% group_by(GENDER,CANCER_TYPE,PATIENT_ID) %>% dplyr::summarise() %>% 
  group_by(GENDER,CANCER_TYPE) %>% dplyr::summarise(N.P = n()) -> patient.numbers.by.cancer

DNA.VAF %>% group_by(GENDER,PATIENT_ID) %>% dplyr::summarise() %>% 
  group_by(GENDER) %>% dplyr::summarise(N.P = n()) -> patient.numbers


p53.mutants = filter(DNA.VAF,HUGO_SYMBOL=='TP53')
p53.mutants = unique(p53.mutants$PATIENT_ID)

filter.out = c("Silent",'Intron','IGR',"In_Frame_Ins" ,"In_Frame_Del", "lincRNA" )

filter(DNA.VAF,!(VARIANT_CLASSIFICATION%in%filter.out)) -> DNA.VAF
filter(DNA.VAF,CHROMOSOME=='X') -> DNA.VAF
DNA.VAF$SUM = DNA.VAF$DNA_ALT_COUNT + DNA.VAF$DNA_REF_COUNT
DNA.VAF %>% mutate(VAF = ifelse(SUM>=10,DNA_ALT_COUNT/SUM,0.0)) %>% filter(SUM>=10) -> DNA.VAF

##Join VAF AND RMAF
RMAF$MUTATION_ID = as.vector(RMAF$MUTATION_ID)
DNA.VAF$MUTATION_ID = as.vector(DNA.VAF$MUTATION_ID)

RMAF %>% left_join(y = DNA.VAF[,c('MUTATION_ID','VAF')],by = c('MUTATION_ID')) -> RMAF


####P53 set
p53.int = read.csv("~/Documents/PhD/GenderAnalysis/TP53_interactions/Candidates.2018.03.01.csv")

p53.high = p53.int[which(p53.int$score>=0.6),'ID']
p53.med = p53.int[which(p53.int$score>=0.4),'ID']
p53.low = p53.int$ID

RMAF$TP53.INT.HIGH = 'No'
RMAF[which(RMAF$HUGO_SYMBOL%in%p53.high),'TP53.INT.HIGH'] = 'Yes'
RMAF$TP53.INT.MED = 'No'
RMAF[which(RMAF$HUGO_SYMBOL%in%p53.med),'TP53.INT.MED'] = 'Yes'
RMAF$TP53.INT.LOW = 'No'
RMAF[which(RMAF$HUGO_SYMBOL%in%p53.low),'TP53.INT.LOW'] = 'Yes'


RMAF$POINT_COL = 'black'
RMAF[which(RMAF$HUGO_SYMBOL%in%p53.low),'POINT_COL'] = 'red'



##Do Plots##

#Compare check the VAF of EM/SM 
DoMutNumViolinPlotGender(data = RMAF[which(RMAF$SM.RMAF == 1),],data.col = "VAF",fill.col = 'GENDER',
                         x.lab = 'GENDER', y.lab = 'VAF',title.txt = 'DNA VAF for SM by Gender', 
                         show.points = T, points.alpha = 0.1)

DoMutNumViolinPlotGender(data = RMAF[which(RMAF$EM.RMAF == 1),],data.col = "VAF",fill.col = 'GENDER', show.points = T,
                         x.lab = 'GENDER', y.lab = 'VAF',title.txt = 'DNA VAF for EM by Gender', points.alpha = 0.1)


DoMutNumViolinPlotGender(data = DNA.VAF,data.col = "VAF",fill.col = 'GENDER',show.points = T,
                         x.lab = 'GENDER', y.lab = 'VAF',title.txt = 'VAF by Gender')

DoMutNumViolinPlotGenderFacet(data = DNA.VAF,data.col = "VAF",fill.col = 'GENDER',
                              facet.col = 'CANCER_TYPE',x.lab = 'GENDER', y.lab = 'VAF',title.txt = 'VAF by Gender')



DoMutNumViolinPlotGender(data = RMAF,data.col = "RMAF",fill.col = 'GENDER', show.points = T, points.alpha = 0.050,
                         x.lab = 'GENDER', y.lab = 'RMAF',title.txt = 'RMAF by Gender' )

DoMutNumViolinPlotGenderFacet(data = RMAF,data.col = "RMAF",fill.col = 'GENDER',
                              facet.col = 'CANCER_TYPE',x.lab = 'GENDER', y.lab = 'RMAF',title.txt = 'RMAF by Gender')


####Log Ratio plots
#SM
sm.table = read.csv(paste(output.dir,'SM.Chi.Test.Propensity.csv',sep = ''))
sm.table = sm.table[which(sm.table$silent.mutations.females+sm.table$silent.mutations.males>=10),]
sm.table$total.silent.mut = sm.table$silent.mutations.females+sm.table$silent.mutations.males
sm.table$log2ratio = 0
sm.table$log2ratio = (sm.table$silent.mutations.females / sm.table$silent.mutations.males) / 
  (sm.table$total.females / sm.table$total.males)

sm.table$log2ratio = ifelse(sm.table$log2ratio == Inf,50,sm.table$log2ratio)
sm.table$log2ratio = log(sm.table$log2ratio)
sm.table$text = ifelse(sm.table$pvalue<0.05,as.vector(sm.table$feature),'')

log2ratiodir = 'log2ratioPlots/'
sm.table$p53.int = 'No'
sm.table[which(sm.table$feature%in%p53.low),'p53.int'] = 'Yes'
sm.table$text.p53 = ifelse(sm.table$pvalue<=0.05&sm.table$p53.int=='Yes',as.vector(sm.table$feature),'')
sm.table$sig = ifelse(sm.table$pvalue<=0.05,'yes','no')
library(ggplot2)
for (d in datasets){

  
  sm.table$p53.int = 'No'
  sm.table[which(sm.table$feature%in%p53.high),'p53.int'] = 'Yes'
  sm.table$feature.show = ifelse(sm.table$pvalue<=0.05&sm.table$p53.int=='Yes',as.vector(sm.table$feature),'')
  
  pos <- position_jitter(width = 0.5,seed = 1)
  p <- ggplot(sm.table[which(sm.table$DATASET==d),],aes(x=log2ratio,y = -log10(pvalue),size=total.silent.mut,colour=p53.int,shape=sig))
  tiff(filename = paste(output.dir,log2ratiodir,d,'.','SM.p53HIGH.log2ratio.pval.NO.LABEL.tiff',sep = ''),width = 2400,height = 1600,res = 200)
  print(
  p + 
    geom_point(position = pos) +
    coord_cartesian(xlim=c(-4,4)) + 
    theme(legend.title = element_text(size=18), legend.text = element_text(size=18)) 
  # +
  #   geom_label_repel(position = pos,aes(fill = p53.int,label=feature.show ),
  #                   size = 2, fontface='bold',color="white",box.padding = 0.35,
  #                    segment.colour = "grey50")
  )
  dev.off()
    
}

#EM
em.table = read.csv(paste(output.dir,'EM.Chi.Test.Propensity.csv',sep = ''))
em.table = em.table[which(em.table$expressed.mutations.females+em.table$expressed.mutations.males>=10),]
em.table$total.expressed.mut = em.table$expressed.mutations.females+em.table$expressed.mutations.males
em.table$log2ratio = 0
em.table$log2ratio = (em.table$expressed.mutations.females / em.table$expressed.mutations.males) / 
  (em.table$total.females / em.table$total.males)

em.table$log2ratio = ifelse(em.table$log2ratio == Inf,30,em.table$log2ratio)
em.table$log2ratio = ifelse(em.table$log2ratio == 0,.05,em.table$log2ratio)

em.table$log2ratio = log(em.table$log2ratio)
em.table$text = ifelse(em.table$pvalue<0.05,as.vector(em.table$feature),'')

log2ratiodir = 'log2ratioPlots/'
em.table$p53.int = 'No'
em.table[which(em.table$feature%in%p53.low),'p53.int'] = 'Yes'
em.table$text.p53 = ifelse(em.table$pvalue<=0.05&em.table$p53.int=='Yes',as.vector(em.table$feature),'')
em.table$sig = ifelse(em.table$pvalue<=0.05,'yes','no')
em.table$feature.show = ifelse(em.table$pvalue<=0.05&em.table$p53.int=='Yes',as.vector(em.table$feature),'')
library(ggplot2)
for (d in datasets){
  
  
  em.table$p53.int = 'No'
  em.table[which(em.table$feature%in%p53.low),'p53.int'] = 'Yes'
  pos <- position_jitter(width = 0.5,seed = 1)
  p <- ggplot(em.table[which(em.table$DATASET==d),],aes(x=log2ratio,y = -log10(pvalue),size=total.expressed.mut,colour=p53.int,shape=sig))
  tiff(filename = paste(output.dir,log2ratiodir,d,'.','EM.p53LOW.log2ratio.pval.LABEL.tiff',sep = ''),width = 2400,height = 1600,res = 200)
  print(
    p + geom_point(position = pos) +
      geom_jitter(position = pos) +
      coord_cartesian(xlim=c(-4,4)) + 
      theme(legend.title = element_text(size=18), legend.text = element_text(size=18)) 
    +
     geom_label_repel(position = pos,aes(log2ratio,-log10(pvalue),fill = p53.int,label=feature.show ),
                     size = 2, fontface='bold',color="white",box.padding = 0.35,
                     segment.colour = "grey50")
  )
  dev.off()
  
}


#DNAmuts
sm.table = read.csv(paste(output.dir,'DNA.muts.Chi.Test.Propensity.csv',sep = ''))
sm.table = sm.table[which(sm.table$mutations.females+sm.table$mutations.males>=10),]
sm.table$total.mut = sm.table$mutations.females+sm.table$mutations.males
sm.table$log2ratio = 0
sm.table$log2ratio = (sm.table$mutations.females / sm.table$mutations.males) / 
  (sm.table$total.females / sm.table$total.males)

sm.table$log2ratio = ifelse(sm.table$log2ratio == Inf,50,sm.table$log2ratio)
sm.table$log2ratio = log(sm.table$log2ratio)
sm.table$text = ifelse(sm.table$pvalue<0.05,as.vector(sm.table$feature),'')

log2ratiodir = 'log2ratioPlots/'
sm.table$p53.int = 'No'
sm.table[which(sm.table$feature%in%p53.high),'p53.int'] = 'Yes'
sm.table$text.p53 = ifelse(sm.table$pvalue<=0.05&sm.table$p53.int=='Yes',as.vector(sm.table$feature),'')
sm.table$sig = ifelse(sm.table$pvalue<=0.05,'yes','no')
library(ggplot2)
for (d in datasets){
  
  
  sm.table$p53.int = 'No'
  sm.table[which(sm.table$feature%in%p53.high),'p53.int'] = 'Yes'
  
  p <- ggplot(sm.table[which(sm.table$DATASET==d),],aes(log2ratio,-log10(pvalue)))
  jpeg(filename = paste(output.dir,log2ratiodir,d,'.','DNA.muts.p53HIGH.log2ratio.pval',sep = ''),width = 6000,height = 4000,res = 400)
  print(
    p + geom_point(aes(size=total.mut,colour=p53.int,shape=sig)) + coord_cartesian(xlim=c(-4,4)) + 
      geom_text(aes(label=feature),size = 2, hjust=0.5, vjust=0) 
  )
  dev.off()
  
}


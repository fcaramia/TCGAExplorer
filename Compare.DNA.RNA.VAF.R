rm(list = ls())
library(dplyr)
library(ggplot2)
library(gtools)
source('ExpressionPlotsFunctions.R')
library(data.table)
#TCGA Expression plots
#Datasets and directories and variables
datasets = c("ESCA","HNSC","LUSC","BLCA","LIHC","STAD","LGG","COAD","PAAD","READ","SKCM","LUAD")
#datasets = c("READ")
output.dir = "~/Documents/PhD/GenderAnalysis/TCGA/Analysis/TCGAExpressionExplorerOutput/"

###Clinical Data####
demo = read.csv("~/Documents/PhD/GenderAnalysis/TCGA/Analysis/all.TCGA.curated.clinical.csv", as.is = T)
demo$PATIENT_ID = toupper(demo$PATIENT_ID)
###############

#READ TP53 Interactors

p53.int = c("EFNB1","PSMD10","STAG2","MAGED2","EDA2R","RNF128","HUWE1","TAF1",
            "BMX","MAGEB18","TSC22D3","MAGEA2B","ACSL4","MAGEH1","PLAC1","BRCC3",
            "SH2D1A","UBE2A","SLC25A5","ATRX","RPS4X","AR","APEX2","FOXP3","RBM3",
            "ARX","POLA1","PHKA2","UTP14A","DDX3X","CUL4B","CITED1","YY2")

###Read DNA VAF#####
DNA.VAF = fread("~/Documents/PhD/GenderAnalysis/TCGA/Analysis/full.reduced.all.raw.TCGA.curated.mutations.csv")
DNA.VAF$PATIENT_ID = gsub("TCGA\\-([[:alnum:]]{2})\\-([[:alnum:]]{4}).*","\\2", DNA.VAF$SAMPLE_ID)
DNA.VAF %>% left_join(demo[,c('GENDER','PATIENT_ID')]) -> DNA.VAF
DNA.VAF %>% group_by(GENDER,CANCER_TYPE,PATIENT_ID) %>% dplyr::summarise() %>% 
  group_by(GENDER,CANCER_TYPE) %>% dplyr::summarise(N.P = n()) -> patient.numbers
p53.mutants = filter(DNA.VAF,HUGO_SYMBOL=='TP53')
p53.mutants = unique(p53.mutants$PATIENT_ID)

##RMAF DATA
RMAF = fread("~/Documents/PhD/GenderAnalysis/TCGA/Analysis/all.TCGA.ExMut.V5.csv")
RMAF$PATIENT_ID = gsub("TCGA\\-([[:alnum:]]{2})\\-([[:alnum:]]{4}).*","\\2", RMAF$SAMPLE)
RMAF$TP53_STATUS = ifelse(RMAF$PATIENT_ID%in%p53.mutants,'Mt','Wt')
RMAF$TP53_INT = ifelse(RMAF$HUGO_SYMBOL%in%p53.int,'I','N')

signature.file =  "~/Documents/PhD/GenderAnalysis/TCGA/Analysis/TP53.lists.2.txt"
###Read Signatures####
sig.lists = read.delim(signature.file, sep = " ")

#Clean low reads
RMAF$SUM = RMAF$REF_COUNT + RMAF$ALT_COUNT
RMAF %>% filter(SUM>10) -> RMAF
RMAF$RMAF = RMAF$ALT_COUNT/RMAF$SUM
RMAF$LOGIT.RMAF = logit(x = RMAF$RMAF,min = 0-0.1 , max = 1+0.1)

#Define Expressed Mutation
cutoff = 0.3
RMAF$EM = ifelse(RMAF$RMAF>=cutoff,1,0)
RMAF$SM = ifelse(RMAF$RMAF>=cutoff,0,1)


RMAF.INT = filter(RMAF,TP53_INT=='I')
RMAF.NO.INT = filter(RMAF,TP53_INT=='N')

source('ExpressionPlotsFunctions.R')
DoDensityPlotFacetGenderRMAF(dat = RMAF, dat.col = 'LOGIT.RMAF', facet.col = 'TP53_INT'
                             ,fill.col = 'GENDER',title.txt = 'RMAF densities by Cancer')

res = data.frame(HUGO_SYMBOL=unique(RMAF$HUGO_SYMBOL),ALL_CANCERS=1)

for (g in unique(RMAF$HUGO_SYMBOL)){
    g.dat = filter(RMAF,HUGO_SYMBOL == g)
    g.dat %>% group_by(GENDER) %>% dplyr::summarise(TEM = sum(EM), TSM = sum(SM)) 
  
    s = tryCatch(
    {
      summary(lm(LOGIT.RMAF~GENDER,data = filter(RMAF,HUGO_SYMBOL==g)))$coeff[1,4] 
    },error = function(e){
      return(1)
    },warning = function(w){
      return(1)
    }
    )
  res[which(res$HUGO_SYMBOL==g),'ALL_CANCERS'] = s
}

for(c in unique(RMAF$CANCER_TYPE)){
  res[,c] = 1
  
  for (g in unique(RMAF$HUGO_SYMBOL)){
    s = tryCatch(
      {
        summary(lm(LOGIT.RMAF~GENDER,data = filter(RMAF,HUGO_SYMBOL==g, CANCER_TYPE==c)))$coeff[1,4] 
      },error = function(e){
        return(1)
      },warning = function(w){
        return(1)
      }
    )
    res[which(res$HUGO_SYMBOL==g),c] = s
  }
  
}

write.csv(res,"~/Documents/PhD/GenderAnalysis/TCGA/Analysis/TCGAExpressionExplorerOutput/AllDataSets/RMAF~Gender.P-values.all.genes.csv",
          row.names = F)

#Do hypergeom tests

for (c in colnames(res)){
  
  if (c=='HUGO_SYMBOL'){
    next
  }
  if (c=='ALL_CANCERS'){
    tot.genes = length(unique(RMAF$HUGO_SYMBOL))
    int.in.tot = p53.int[p53.int%in%unique(RMAF$HUGO_SYMBOL)]
    tot.int = length(int.in.tot)
    
  }else{
    tot.genes = length(unique(filter(RMAF,CANCER_TYPE==c)$HUGO_SYMBOL))
    int.in.tot = p53.int[p53.int%in%unique(filter(RMAF,CANCER_TYPE==c)$HUGO_SYMBOL)]
    tot.int = length(int.in.tot)
  }
  
  
  sig.genes = as.vector(res[which(res[,c]<=0.05),'HUGO_SYMBOL'])
  sig.int = length(p53.int[p53.int%in%sig.genes])
  
  p.val = phyper(sig.int,tot.int,tot.genes-tot.int,length(sig.genes),lower.tail = F)
  
  print(c)
  print(paste('Total Genes ',tot.genes))
  print(paste('Interactors present ', tot.int))
  print(paste(int.in.tot,collapse = ' '))
  print(paste('Sig Interactors ', sig.int))
  print(paste(p53.int[p53.int%in%sig.genes],collapse = ' '))
  print(paste('Sig Genes ', length(sig.genes)))
  print(paste('p-val',p.val))
  
}


###Read DNA VAF#####
DNA.VAF = fread("~/Documents/PhD/GenderAnalysis/TCGA/Analysis/full.reduced.all.raw.TCGA.curated.mutations.csv")

###Read RMAF
RMAF = fread("~/Documents/PhD/GenderAnalysis/TCGA/Analysis/TCGAExpressionExplorerOutput/AllDataSets/Computed.All.RMAF.csv")

RMAF$SUM = RMAF$REF_COUNT + RMAF$ALT_COUNT

RMAF %>% mutate(RMAF = ifelse(SUM>=1,ALT_COUNT/SUM,0.0)) %>% filter(!(EXPRS.GENE=='NOT.EXPRS'&SUM<10)) %>%
  filter(SUM>10) %>%  filter(CANCER_TYPE%in%datasets)-> RMAF


#Clean low reads
DNA.VAF$SUM = DNA.VAF$DNA_ALT_COUNT + DNA.VAF$DNA_REF_COUNT
DNA.VAF %>% filter(SUM>10) -> DNA.VAF
DNA.VAF$DNA.VAF = DNA.VAF$DNA_ALT_COUNT/DNA.VAF$SUM
DNA.VAF$LOGIT.DNA.VAF = logit(x = DNA.VAF$DNA.VAF,min = 0-0.1 , max = 1+0.1)
DNA.VAF$PATIENT_ID = gsub("TCGA\\-([[:alnum:]]{2})\\-([[:alnum:]]{4}).*","\\2", DNA.VAF$SAMPLE_ID)
DNA.VAF$MUTATION_ID = paste(DNA.VAF$PATIENT_ID,DNA.VAF$HUGO_SYMBOL,DNA.VAF$CHROMOSOME,
                            DNA.VAF$START_POSITION,DNA.VAF$END_POSITION,DNA.VAF$TUMOR_SEQ_ALLELE2,sep = '.')

DNA.VAF %>% filter(MUTATION_ID%in%RMAF$MUTATION_ID) %>% 
  select(MUTATION_ID,DNA.VAF) -> DNA.comp

RMAF %>% select(MUTATION_ID,RMAF,GENDER) -> RNA.comp

  
###Merge DNA VAF and RMAF 

DNA.comp %>% left_join(RNA.comp) -> COMB

comb.test = COMB[COMB$GENDER=='female',]
cor.test(comb.test$DNA.VAF,comb.test$RMAF)



COMB$EFFECT = ifelse(COMB$DNA.VAF>0.3&COMB$RMAF>0.3,'Normal',
                     ifelse(COMB$DNA.VAF>0.3&COMB$RMAF<0.3,'Silenced',
                            ifelse(COMB$DNA.VAF<0.15&COMB$RMAF<0.15,'Noise','Clonality/LOH')))

COMB %>% filter(CHROMOSOME!='X'&CHROMOSOME!='Y') -> COMB.AUTO
COMB %>% filter(CHROMOSOME=='X') -> COMB.X

#COMB %>% select(GENDER,LOGIT.DNA.VAF,LOGIT.RMAF) %>% melt() -> MELTED

#Scattered RMAF vs DNA VAF

p <- ggplot(COMB.X,aes(RMAF,DNA.VAF)) 
p <- p +  geom_point(size = 0.2, alpha = 1/20) + facet_grid(CANCER_TYPE~GENDER) + coord_fixed()
p + stat_density_2d(aes(fill = ..level..), geom="polygon")+
  scale_fill_gradient(low="blue", high="red")


p <- ggplot(COMB.AUTO,aes(RMAF,DNA.VAF)) 
p <- p +  geom_point(size = 0.2, alpha = 1/20) + facet_grid(CANCER_TYPE~GENDER) + coord_fixed()
p + stat_density_2d(aes(fill = ..level..), geom="polygon")+
  scale_fill_gradient(low="blue", high="red")



#GEnome wide RMAF
COMB$EFFECT = ifelse(COMB$DNA.VAF>0.3&COMB$RMAF>0.3,'Normal',
                     ifelse(COMB$DNA.VAF>0.3&COMB$RMAF<0.3,'Silenced',
                            ifelse(COMB$DNA.VAF<0.15&COMB$RMAF<0.15,'Normal','Normal')))
COMB %>% filter(CHROM!='X'&CHROM!='Y') -> COMB.AUTO
COMB %>% filter(CHROM=='X') -> COMB.X

p <- ggplot(RMAF[which(RMAF$CHROM=='X'),],aes(START_POSITION,RMAF)) 
p +  geom_point(size = 1.5, alpha=0.05) + facet_grid(CANCER_TYPE~GENDER,scales =  'free_x') 


#MELTED$VAR = paste(MELTED$GENDER,MELTED$variable,sep = '.') 

DoDensityPlotGenderRMAF(dat = MELTED, dat.col = 'value' ,
                        fill.col = 'VAR', 
                        title.txt = paste('LIHC RMAF and DNA VAF'))




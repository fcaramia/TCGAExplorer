library(data.table)
pats = fread('~/Documents/PhD/GenderAnalysis/TCGA/Analysis/Normal.Samples.in.SBG.csv')
germs = fread('~/Documents/PhD/GenderAnalysis/TCGA/Analysis/Germline.Normals.ExMut.csv')
pats$PATIENT_ID = gsub("TCGA\\-([[:alnum:]]{2})\\-([[:alnum:]]{4}).*","\\2", pats$sample_id)
pats$SAMPLE_ID = pats$sample_id                       
pats =  pats[,c('SAMPLE_ID','PATIENT_ID')]

write.csv(pats, '~/Documents/PhD/GenderAnalysis/TCGA/Analysis/Normal.Patients.in.SBG.csv', row.names = F )

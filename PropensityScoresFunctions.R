source('PropensityScores/cal.R')

DoCalcWeights <- function(clinical.data, default.confounding.factors, extra.cofounding.factors)
{
  cl = as.data.frame(clinical.data)
  cfe = extra.cofounding.factors
  cfd = default.confounding.factors
  nr.patients = dim(cl)[1]
  if(nr.patients<=60)
  {
    return(NULL)  
  }
    
  analysis = 'gender'
  badchars <- "[\xb5]|[\n]|[,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[/]|[\\]|[.]|[_]|[ ]"
  
  for (i in unique(c(cfd,cfe)))
  {
    
    cutoff = ceiling(nr.patients*.1)
    nas = length(cl[which(is.na(cl[,i])),i])
    if(nas>cutoff){
      cfe = cfe[cfe!= i]
    }
    else{
      cl=cl[which(!is.na(cl[,i])),]
      
    }
  }
  
  #Calculate samples needed to use propensity scores.
  cfd.domain = 20
  cfe.domain = 0
  for(c in cfe){
    cfe.domain = cfe.domain + length(levels(as.factor(cl[,c])))
  }
  cfe.domain = cfe.domain*10
  if (nr.patients <= cfe.domain)
  {
    cfe = NULL
  }
  cf = unique(c(cfd,cfe))
  
  cl = cl[,c(cf)]
  cl$GENDER <- ifelse(cl$GENDER=="female",1,0)
  colnames(cl)[which(colnames(cl)=='GENDER')] <- "Z"
  rownames(cl) = cl$PATIENT_ID
  cl$PATIENT_ID = NULL
  
  #Calculate Scores for all samples 
  # convert clinical variables to dummy
  library(dummies)
  dummy.feature <- setdiff(colnames(cl),c("Z",'AGE_AT_INITIAL_PATHOLOGIC_DIAGNOSIS'))
  if (!is.null(cfe))
  {
    cl <- dummy.data.frame(cl, names=cfe)
  }
  dummy.list <- attr(cl,"dummies")
  rm.col <- c()
  for (i in 1:length(dummy.list))
  {
    rm.col <- c(rm.col, dummy.list[[i]][length(dummy.list[[i]])])
  }
  if (!is.null(rm.col))
  {
    cl <- cl[,-rm.col]
  }
  
  cl$X0 <- rep(1, nrow(cl))
  exclude.col <- match(c("Z","X0"), colnames(cl))
  colnames(cl) <- toupper(gsub(pattern=badchars, replacement=".", x=colnames(cl)))
  #colnames(cl) <- gsub(" ", ".", colnames(cl))
  form <- as.formula(paste("Z~",paste(colnames(cl)[-exclude.col],collapse="+"),sep=""))
  ##### End of processing clinical data
  
  ############################
  ##### perform calculation
  
  #Mutation.pri will be the input for the analysis. Each line is a sample while each column is a gene
  wt.all <- weight.test(cl, form, NULL, just.calc.weights = T, is.continuous=TRUE,weight=ifelse(analysis=="gender","MW","ATT"),mirror.plot=FALSE, cancer, "Dummy", outdir=NULL)
  cl$WEIGHTS = wt.all
  cl$PATIENT_ID = rownames(cl)
  return(cl)
}

DoStatTest <- function(clinical.data, molecular.data, default.confounding.factors, extra.cofounding.factors, is.continuous=F){
  cl = as.data.frame(clinical.data)
  cfe = extra.cofounding.factors
  cfd = default.confounding.factors
  nr.patients = dim(cl)[1]
  if(nr.patients<=60)
  {
    return(NULL)  
  }
  
  analysis = 'gender'
  badchars <- "[\xb5]|[\n]|[,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[/]|[\\]|[.]|[_]|[ ]"
  
  
  for (i in unique(c(cfd,cfe)))
  {
    
    cutoff = ceiling(nr.patients*.1)
    nas = length(cl[which(is.na(cl[,i])),i])
    if(nas>cutoff){
      cfe = cfe[cfe!= i]
    }
    else{
      cl=cl[which(!is.na(cl[,i])),]
      
    }
  }
  
  #Calculate samples needed to use propensity scores.
  cfd.domain = 20
  cfe.domain = 0
  for(c in cfe){
    cfe.domain = cfe.domain + length(levels(as.factor(cl[,c])))
  }
  cfe.domain = cfe.domain*10
  if (nr.patients <= cfe.domain)
  {
    cfe = NULL
  }
  if(length(cfe)==0)
  {
    cfe=NULL
  }
  cf = unique(c(cfd,cfe))
  
  cl = cl[,c(cf)]
  cl$GENDER <- ifelse(cl$GENDER=="female",1,0)
  colnames(cl)[which(colnames(cl)=='GENDER')] <- "Z"
  rownames(cl) = cl$PATIENT_ID
  cl$PATIENT_ID = NULL
  
  #Calculate Scores for all samples 
  # convert clinical variables to dummy
  library(dummies)
  dummy.feature <- setdiff(colnames(cl),c("Z",'AGE_AT_INITIAL_PATHOLOGIC_DIAGNOSIS'))
  if (!is.null(cfe))
  {
    cl <- dummy.data.frame(cl, names=cfe)
  }
  dummy.list <- attr(cl,"dummies")
  rm.col <- c()
  for (i in 1:length(dummy.list))
  {
    rm.col <- c(rm.col, dummy.list[[i]][length(dummy.list[[i]])])
  }
  if (!is.null(rm.col))
  {
    cl <- cl[,-rm.col]
  }
  
  cl$X0 <- rep(1, nrow(cl))
  exclude.col <- match(c("Z","X0"), colnames(cl))
  colnames(cl) <- toupper(gsub(pattern=badchars, replacement=".", x=colnames(cl)))
  #colnames(cl) <- gsub(" ", ".", colnames(cl))
  form <- as.formula(paste("Z~",paste(colnames(cl)[-exclude.col],collapse="+"),sep=""))
  ##### End of processing clinical data
  
  ############################
  ##### perform calculation
  
  #Mutation.pri will be the input for the analysis. Each line is a sample while each column is a gene
  wt.all <- weight.test(data = cl, form=form, molecular.pri = molecular.data,just.calc.weights = F, is.continuous=is.continuous,weight=ifelse(analysis=="gender","MW","ATT"),mirror.plot=FALSE, cancer, "Dummy", outdir=NULL)
  #cl$WEIGHTS = wt.all
  #cl$PATIENT_ID = rownames(cl)
  return(wt.all)
  
}
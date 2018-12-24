# Yuan Yuan, 2015-03-11 weighted test agains continuous response variables, such as gene expression
# Hu Chen, maintaining and updating

library(ggplot2)
library(grid)
library(gridExtra)
library(gtable)
library(weights)

hist.mirror2 <- function(data.plot, label)
    {
        myplot <- ggplot(data.plot, aes(ps, fill=status))+geom_histogram(data=subset(data.plot, status=="Yes"), aes(ps, y=..count..))+
            geom_histogram(data=subset(data.plot,status=="No"), aes(ps, y= -..count..))+
                scale_fill_hue(label)+
                    labs(x="Propensity Score")
        print(myplot)
    }


weight.test <- function(data,form, molecular.pri, just.calc.weights = T, is.continuous=TRUE, weight="ATT", mirror.plot=FALSE,cancer, data.type, outdir=".", perm=FALSE, seed=seed)
    {        
        common <- intersect(rownames(data), rownames(molecular.pri))
        print(paste("Number of samples:", dim(data)[1]))
        
        if (!just.calc.weights){
          molecular.common <- molecular.pri[match(common, rownames(molecular.pri)),]
          clinical.common <- data[match(common, rownames(data)),]
          
        }
        else
        {
          clinical.common <- data
        }
        #clinical.common <- data
        print(summary(factor(clinical.common$Z)))
        # write sample info to file
        #write(rbind(rownames(clinical.common),ifelse(clinical.common$Z==1, "Female","Male")), file=paste("sample_list/",cancer,"_",data.type,"_samples.txt", sep=""),ncol=2)
        
        print("----------------")

        age.cutoff <- 65
        print(paste("age <", age.cutoff,":", length(which(clinical.common$AGE.AT.INITIAL.PATHOLOGIC.DIAGNOSIS< age.cutoff))))
        print(paste("age >=", age.cutoff,":", length(which(clinical.common$AGE.AT.INITIAL.PATHOLOGIC.DIAGNOSIS>= age.cutoff))))
        print("----------------")
        

        if(perm)
            {
                n <- nrow(clinical.common)
                set.seed(seed)
                perm <- sample(1:n,n)
                clinical.common$Z <- clinical.common$Z[perm]
            }
 
        print(paste("Weighting scheme:", weight))
        source("PropensityScores/GIPW_function_omega.R")
        ans <- GIPW.std.omega(dat=clinical.common, form.ps=form, weight=weight,trt.est=F)
        # draw mirror plot
        
        if(mirror.plot)
            {
                ps <- ans$ps
                status <- ifelse(clinical.common$Z==1, "Yes","No")
                data.plot <- data.frame(ps, status )
                pdf(paste(outdir,"/",cancer,"_", data.type, ifelse(perm,paste("_perm_",seed, sep=""), ""),"_mirror_raw.pdf",sep=""))
                if(analysis=="omt")
                    {
                        label="OMT"
                    }
                if(analysis=="gender")
                    {
                        label="female"
                    }
                if(analysis=="race")
                    {
                        label="White"
                    }
                hist.mirror2(data.plot, label)
                dev.off()
            }
        wt <- ans$W
        source("PropensityScores/check_balance.R")
        index.tr <- which(clinical.common$Z==1)
        index.ctr <- which(clinical.common$Z==0)

        if(mirror.plot)
        {
          ps <- ans$ps
          ps <- ps*wt
          status <- ifelse(clinical.common$Z==1, "Yes","No")
          data.plot <- data.frame(ps, status )
          pdf(paste(outdir,"/",cancer,"_", data.type, ifelse(perm,paste("_perm_",seed, sep=""), ""),"_mirror_WT.pdf",sep=""))
          if(analysis=="omt")
          {
            label="OMT"
          }
          if(analysis=="gender")
          {
            label="female"
          }
          if(analysis=="race")
          {
            label="White"
          }
          hist.mirror2(data.plot, label)
          dev.off()
        }
        
        cutoff=0.1 # 10%
        for (i in 2:(ncol(clinical.common)-1))
            {
                print(paste("check", colnames(clinical.common)[i]))
                #print(summary(clinical.common[,i]))
                std.diff <- check.balance(index.tr, index.ctr, clinical.common[,i], wt, printout=TRUE)
                if(std.diff>cutoff)
                    {
                        print(paste("Fail: standardized difference= ", std.diff))
                        print("=======================================")
                        stop()
                    }else{
                        print("Pass!")
                    }
            }
        
        if (just.calc.weights){
          return (wt)  
        }
        
        molecular.pvalues <- c()
        molecular.coefs <- c()
        molecular.dfs <- c()
        molecular.chisqs <- c()
        molecular.0 <- c()
        molecular.1 <- c()
        molecular.0.w <- c()
        molecular.1.w <- c()
        for (i in 1:ncol(molecular.common))
            {
                tmp <- by(molecular.common[,i], clinical.common$Z, mean, na.rm = T)
                molecular.0 <- c(molecular.0, tmp[1])
                molecular.1 <- c(molecular.1, tmp[2])
                # get weighted mean
                molecular.0.w <- c(molecular.0.w, sum(molecular.common[index.ctr,i]*wt[index.ctr])/sum(wt[index.ctr]))
                molecular.1.w <- c(molecular.1.w, sum(molecular.common[index.tr,i]*wt[index.tr])/sum(wt[index.tr]))
                #print(i)
                if (i%%1000==0)
                    {
                        print(i)
                    }
                
                if(is.continuous)
                    {
                        #pvalue <- try(summary(lm(clinical.common$Z~molecular.common[,i], weights=wt))$coef[2,4])
                        pvalue <- try(summary(lm(molecular.common[,i]~clinical.common$Z, weights=wt))$coef[2,4])
                        coef <- try(summary(lm(molecular.common[,i]~clinical.common$Z, weights=wt))$coef[2,1])
                        chisq = NA
                        if (class(pvalue)=="try-error")
                            {
                                pvalue <- NA
                                coef <- NA
                            }
                    }else{
                        #Using logistic regression
                        #pvalue <- summary(glm(clinical.common$Z~molecular.common[,i], family=binomial, weights=wt))$coef[2,4]
                        #pvalue <- summary(glm(molecular.common[,i]~clinical.common$Z, family=binomial, weights=wt))$coef[2,4]
                        #coef <- summary(glm(molecular.common[,i]~clinical.common$Z, family=binomial, weights=wt))$coef[2,1]

						            #Using weighted chi-square test
					             	wtd.test <- wtd.chi.sq(molecular.common[,i],clinical.common$Z,weight=wt)
                        pvalue <- as.numeric(wtd.test["p.value"])
					            	wcst <- by(molecular.common[,i]*wt, clinical.common$Z, mean)
						            df <- as.numeric(wtd.test["df"])
						            chisq <- as.numeric(wtd.test["Chisq"])
					            	coef <- wcst[2] - wcst[1]
                        if (class(pvalue)=="try-error")
                            {
                                pvalue <- NA
                                coef <- NA
                                chisq <- NA
                                df <- NA
                            }
                    }
                    if(is.na(molecular.0[i])|is.na(molecular.1[i])){
                      molecular.pvalues <- c(molecular.pvalues, pvalue)
                    }else{
                      if( (molecular.0[i] == 0 && molecular.1[i] == 0))
                      {
                        molecular.pvalues <- c(molecular.pvalues,1)
                      }else{
                        molecular.pvalues <- c(molecular.pvalues, pvalue)
                      }  
                    }
		                
                    molecular.dfs <- c(molecular.dfs, df)
                    molecular.chisqs <- c(molecular.chisqs, chisq)
                    molecular.coefs <- c(molecular.coefs, coef)
            }
            molecular.fdr <- p.adjust(molecular.pvalues,"BH")
            if(is.continuous==F){
              return(list(feature=colnames(molecular.common),
                          chisq=molecular.chisqs,df=molecular.dfs,
                          coef=molecular.coefs, pvalue=molecular.pvalues,
                          fdr=molecular.fdr , mean.0= molecular.0, mean.1= molecular.1, 
                          mean.0.w=molecular.0.w, mean.1.w=molecular.1.w))  
            }
            else
            {
              return(list(feature=colnames(molecular.common),
                          coef=molecular.coefs, pvalue=molecular.pvalues,
                          fdr=molecular.fdr , mean.0= molecular.0, mean.1= molecular.1, 
                          mean.0.w=molecular.0.w, mean.1.w=molecular.1.w)) 
            }
            
    }

summarize.fdr <- function(molecular.data, molecular.result, cutoff=0.05, print=FALSE )
    {
        pvalue <- molecular.result$pvalue
        fdr <- molecular.result$fdr
        coef <- molecular.result$coef
        mean.0 <- molecular.result$mean.0
        mean.1 <- molecular.result$mean.1
        mean.0.w <- molecular.result$mean.0.w
        mean.1.w <- molecular.result$mean.1.w
        signif.index <- which(fdr<cutoff)
        print(paste("Features with FDR <", cutoff, "=", length(signif.index)))
        if(print)
            {
                print(cbind(colnames(molecular.data)[signif.index],signif(pvalue[signif.index],3), signif(fdr[signif.index],3), signif(coef[signif.index],3), signif(mean.0[signif.index],3), signif(mean.1[signif.index],3)))
            }
        return(list(feature.sig=colnames(molecular.data)[signif.index], pvalue.sig=pvalue[signif.index], fdr.sig= fdr[signif.index], n.sig=length(signif.index), coef.sig=coef[signif.index], mean0.sig=mean.0[signif.index], mean1.sig=mean.1[signif.index], mean0.sig.w=mean.0.w[signif.index], mean1.sig.w=mean.1.w[signif.index]))#, n.sig=length(which(molecular.result$pvalue<cutoff))) )
        
    }


write.summary <- function(summary, cancer, analysis, type)
    {
        file <- paste(cancer,"_",analysis,"_summary_",type,".txt", sep="")
        data <- data.frame(summary$feature.sig, summary$fdr.sig, summary$coef.sig, summary$mean0.sig, summary$mean1.sig, summary$mean0.sig.w, summary$mean1.sig.w)
        if (analysis=="gender")
            {
                colnames(data) <- c("feature","fdr","coef","mean_MALE", "mean_FEMALE", "mean_MALE_weighted", "mean_FEMALE_weighted")
            }
        if(analysis=="race")
            {
                colnames(data) <- c("feature","fdr","coef","mean_nonWHITE", "mean_WHITE","mean_nonWHITE_weighted", "mean_WHITE_weighted")
            }
        if (analysis=="omt")
            {
                colnames(data) <- c("feature","fdr","coef","mean_nonOMT", "mean_OMT","mean_nonOMT_weighted", "mean_OMT_weighted")
            }
        write.table(data,file=file, quote=FALSE, sep="\t", col.names=TRUE, row.names=FALSE)
    }

write.result <- function(result, cancer, analysis, type, outdir,fdr.cutoff=0.05)
    {
        file <- paste(cancer,"_",analysis,"_result_",type,".txt", sep="")
        data <- data.frame(result$feature, result$pvalue, result$fdr, result$coef, result$mean.0, result$mean.1, result$mean.0.w, result$mean.1.w)
        if (analysis=="gender")
            {
                colnames(data) <- c("feature","pvalue","fdr","coef","mean_MALE", "mean_FEMALE", "mean_MALE_weighted", "mean_FEMALE_weighted")
            }
        if(analysis=="race")
            {
                colnames(data) <- c("feature","pvaue","fdr","coef","mean_nonWHITE", "mean_WHITE","mean_nonWHITE_weighted", "mean_WHITE_weighted")
            }
        if (analysis=="omt")
            {
                colnames(data) <- c("feature","fdr","coef","mean_nonOMT", "mean_OMT","mean_nonOMT_weighted", "mean_OMT_weighted")
            }
        write.table(data,file=paste(outdir,"/",file,sep=""), quote=FALSE, sep="\t", col.names=TRUE, row.names=FALSE)      
    }


plot.perm <- function(perm, n, cancer, analysis, type, cutoff)
    {        
        p=length(which(perm>= n))/length(perm)
        print(paste("Permutation P-value =", p))
        print(paste(type,": n.sig =", n))
        print(paste("Median perm n.sig =", median(perm)))
        file <- paste(cancer,"_",analysis,"_perm_",type,"_",cutoff,".pdf", sep="")
        pdf(file, width=3, height=3)
        par(mgp=c(2,1,0))
        if(type=="mut") {main="Mutation"}
        if(type=="cnv") {main="SCNA"}
        if(type=="methy") {main="Methy"}
        if(type=="mRNAseq") {main="mRNA"}
        if(type=="miRNA") {main="miRNA"}
        if(type=="rppa") {main="Protein"}
        if (n > max(perm))
            {
                myhist <- hist(perm, xlim=c(0, max(max(perm),n)), main="", xlab="# features" )
                
            }else{
                myhist <- hist(perm, main="", xlab="# features")
            }
        text(n, max(myhist$counts),paste("p-value =", p), pos=ifelse(n> 0.5* max(perm),2,4), col=ifelse(p<=0.05, "red","blue"), cex=1.1, xpd=T)
        abline(v=n, col=ifelse(p<=0.05, "red","blue"), lty=2, lwd=2)
        dev.off()
    }

perm.cal <- function(cancer, analysis, type, molecular.pri, cutoff=0.05, seedV=1:00)
    {
        print(paste("Calculate permutation for", cancer,":",type,":"))
        load(file=paste(scripts.dir,"/",cancer,"_",analysis,"/",cancer,"_", analysis,"_result.RData",sep=""))
        result <- get(paste(type,".result", sep=""))
        summary <- summarize.fdr(molecular.pri, result, cutoff=cutoff )
        n <- summary$n.sig
        perm <- c()
        print("For permutation:")
        for (seed in seedV)
            {
                ##print(seed)
                perm.result <- perm.summary <- c()
                 if(type=="mature")
                     {
                         load(paste(scripts.dir, "/",cancer,"_",analysis,"/perm_sig_result_mature_",seed,".RData", sep=""))
                     }else if(type=="pre")
                          {
                              load(paste(scripts.dir, "/",cancer,"_",analysis,"/perm_sig_result_pre_",seed,".RData", sep=""))
                          }else{
                
                              load(paste(scripts.dir, "/",cancer,"_",analysis,"/perm_sig_result_",seed,".RData", sep=""))
                          }
                perm.result <- get(paste("perm.",type,".result", sep=""))
                perm.summary <- summarize.fdr(molecular.pri, perm.result, cutoff=cutoff )
                perm <- c(perm, perm.summary$n.sig)
                do.call(rm, list(paste("perm.",type,".result", sep="")))
                               
            }
        ##print (perm)
        ##print (n)
        plot.perm(perm, n, cancer, analysis, type, cutoff)
    }
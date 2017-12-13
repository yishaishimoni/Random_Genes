# Licensed Materials - Property of IBM 5748-XX8
# (C) Copyright IBM Corp. 1992, 1993 All Rights Reserved
# US Government Users Restricted Rights - Use, duplication or disclosure restricted by GSA ADP Schedule Contract with IBM Corp.

# use random sets to predict survival

singleRandomSet <- function(dat, surv, N=50, toPlot=F, subtypes=2, q=0.5, rand=FALSE, replace=TRUE, genes=NULL, 
                            returnGenes=FALSE){
  # if replace=FALSE and N is the number of genes then all the genes are chosen
  tryCatch({
    if(rand){
      pred <- runif(1:nrow(dat))
    } else {
      if(is.null(genes)){
        genes <- sample(x = 1:ncol(dat), size = N, replace = replace)
      }
      tempMat <- dat[ , genes]
      pc <- prcomp(tempMat)
      pred <- pc$x[ ,1]
    }
    # if q is not 0.5 then order the groups according to their average survival
    group <- pred>quantile(pred, probs = q)
    group2 <- pred>quantile(pred, probs = 1-q)
    if (wilcox.test(surv[,1]~group)$statistic > wilcox.test(surv[,1]~group2)$statistic)
      formula <- surv~group
    else
      formula <- surv~group2
    sf <- survfit(formula)

    pVal <- pchisq(q=survdiff(formula)$chisq,df=1,lower.tail=F)
  },
  error= function(e) {
    stop(paste(pred, collapse = ', '))
  })
  if(toPlot){ # this is used for debugging and for individual plots
    plot(sf,lwd=2,col=c("red","blue"),xlab="Time (days)",ylab="% survival",
         main=paste0("Pval=",signif(pVal,digits=3)))
    legend(x="topright",bty="n",legend=names(sf$strata),lwd=2,col=c("red","blue"))
  }
  if(returnGenes)
    return(list(pval=pVal, genes=genes))
  return (pVal)
}

multipleRandomSets <- function(dataObj, d, reps=5000, N=50, q=0.5){
  pVals <- sapply(1:reps, FUN = function(x) singleRandomSet(dataObj$merged.dat, dataObj$Surv, 
                                                            N=N, q=q))
  nullpVals <- sapply(1:reps, FUN = function(x) singleRandomSet(dataObj$merged.dat, dataObj$Surv, 
                                                                N=N, q=q, rand=TRUE))
  
  hist(pVals,30,main=paste("Histogram for pVals for",d))
  
  qqplot(x=nullpVals ,y=pVals, type="l", lwd=2, xlab="null quantiles",
         main=paste("quantile plot for",d))
  abline(a=0, b=1, col="red", lwd=2)
  abline(h=0.05, lwd=2, col="gray")
  
  quants <- exp(seq(log(1/100),log(1),length=100))
  plot(x=quants,y=quantile(pVals, quants), type="l", lwd=2, log="xy", xlab="theoretical quantiles",
       main=paste("quantile plot for",d,"\n(log-scale)"), ylab="pVals")
  abline(a=0, b=1, col="red", lwd=2)
  abline(h=0.05, lwd=2, col="gray")
  
  plot(ecdf(pVals), pch='.', xlab='P-value', ylab='cumulative probability',
       main=paste('Cumulative p-value distibution for',d))
  plot(ecdf(nullpVals), pch='.', col='grey', add=T)
  legend(x='bottomright', legend = c('Random gene-sets','Null distribution'), 
         col=c('black','grey'), lty=1 )
  
  return (list(pvals=pVals, null=nullpVals))
}

random_set_within_cluster <- function(disease, ratios, C, dataObj){
  source('remove_meta_PCNA.r')
  # disease is a string of the TCGA abbreviation. ratios is a list that will be populated with the ratios.
  # C is a list in which C$C is a numeric assignment of samples in dataObj$merged.data into clusters
  for (remove_PCNA in c(FALSE,TRUE)){
    if (remove_PCNA){
      disease <- paste0(disease,'noPCNA')
    }
    ratios[[disease]] <- list()
    for (cluster in 1:max(C$C)){
      ratios[[disease]][[cluster]] <- NaN
      I <- C$C==cluster
      message(paste('checking cluster', cluster, 'using', sum(I), 'samples'))
      if(sum(I) < 20){
        message('skipping - insufficient samples')
        next
      }
      dataObj_C <- dataObj
      dataObj_C$clinical <- dataObj_C$clinical[I]
      dataObj_C$merged.dat <- dataObj_C$merged.dat[I, ]
      dataObj_C$Surv <- dataObj_C$Surv[I, ]
      if (remove_PCNA){
        dataObj_C$merged.dat <- remove_meta_PCNA(dataObj_C)
      }
      pvals <- multipleRandomSets(dataObj=dataObj_C, 
                                  paste0(disease, ' (subclass ', cluster,')'))
      ratios[[disease]][[cluster]] <- sum(pvals$pvals<0.05)/length(pvals$pvals)
      message(ratios[[disease]][[cluster]])
    }
  }
  return(ratios)
}
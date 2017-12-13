# Licensed Materials - Property of IBM 5748-XX8
# (C) Copyright IBM Corp. 1992, 1993 All Rights Reserved
# US Government Users Restricted Rights - Use, duplication or disclosure restricted by GSA ADP Schedule Contract with IBM Corp.

singleSurvivalPrediction <- function(expmat,survObj,method=c("rPCA","random","PCA")[1],scale=F,plot=F,n=50){
  if(method=="rPCA"){
    # choose a random set of genes and use it to predict survival
    rGenes <- sample(1:nrow(expmat),size=n,replace=F)
    if(ncol(expmat)>n){
      pc <- princomp(t(expmat[rGenes,]))
      pred <- pc$scores[,1]
    } else {
      pc <- princomp(expmat[rGenes,])
      pred <- pc$loadings[,1]
    }
  }
  if(method=="random"){
    pred <- runif(n=ncol(expmat))
  }
  if(method=="PCA"){
    pc <- princomp(expmat)
    pred <- pc$loadings[,1]
  }
  formula <- clinical$surv~pred>median(pred)
  sf <- survfit(formula) 
  pVal <- pchisq(q=survdiff(formula)$chisq,df=1,lower.tail=F)
  if(plot){
    plot(sf,lwd=2,col=c("red","blue"),xlab="Time (days)",ylab="% survival",
         main=paste0("Pval=",signif(pVal,digits=3)))
    legend(x="topright",bty="n",legend=names(sf$strata),lwd=2,col=c("red","blue"))
  }
  return(pVal)
}

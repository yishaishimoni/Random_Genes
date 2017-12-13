# Licensed Materials - Property of IBM 5748-XX8
# (C) Copyright IBM Corp. 1992, 1993 All Rights Reserved
# US Government Users Restricted Rights - Use, duplication or disclosure restricted by GSA ADP Schedule Contract with IBM Corp.

# clean up and create a survival object
require(survival)
cleanData <- function(dataObj){
  # take care of objects with multiple data sets
  if (class(dataObj$merged.dat) == "list"){
    merged.dat <- dataObj$merged.dat[[1]]
    for (i in 2:length(dataObj$merged.dat)){
      temp <- dataObj$merged.dat[[i]]
      extraSamples <- setdiff(rownames(temp), rownames(merged.dat))
      merged.dat <- rbind(merged.dat, temp[extraSamples, ])
    }
    dataObj$merged.dat <- merged.dat
  }
  
  dataObj$Surv <- Surv(time = dataObj$merged.dat[,"OS"], event = dataObj$merged.dat[,"status"])
  rownames(dataObj$merged.dat) <- dataObj$merged.dat[ ,1]
  dataObj$merged.dat <- dataObj$merged.dat[ ,-c(1:3)]
  
  # countData <- as.matrix(dataObj$merged.dat)
  # condition <- factor(rep("a",ncol(countData)))
  # dds <- DESeqDataSetFromMatrix(countData, DataFrame(condition), ~condition)
  
  informative <- apply(X = dataObj$merged.dat, MARGIN = 2, FUN = function(x){
    if (any(is.na(x)))
      return (0)
    if (max(table(x))/length(x) > 0.7)
      return (0)
    return (sd(x)/mean(x))
  }) # this will remove genes with NAs, genes that mostly have the same value, and genes with no variance
  
  
  dataObj$merged.dat <- dataObj$merged.dat[ ,informative>0]
  dataObj$merged.dat <- log(1+dataObj$merged.dat)
  return(dataObj)
}
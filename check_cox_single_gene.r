# Licensed Materials - Property of IBM 5748-XX8
# (C) Copyright IBM Corp. 1992, 1993 All Rights Reserved
# US Government Users Restricted Rights - Use, duplication or disclosure restricted by GSA ADP Schedule Contract with IBM Corp.

require(survival)
source('cleanData.r')

check_cox <- function(dataset){
  print(dataset)
  load(dataset)
  dataObj <- cleanData(dataObj)
  pvals <- sapply(X = dataObj$merged.dat[, 1:100], FUN = function(x) summary(coxph(dataObj$Surv~x))$logtest[3])
  return (sum(pvals < 0.050) / length(pvals))
}

diseases <- c("ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "COADREAD", 
              "DLBC", "ESCA", "GBM", "GBMLGG", "HNSC", "KICH", "KIPAN", "KIRC",
              "KIRP", "LAML", "LGG", "LIHC", "LUAD", "LUSC", "MESO", "OV", 
              "PAAD", "PCPG", "PRAD", "READ", "SKCM", "TGCT", "THCA", 
              "THYM", "UCEC", "UCS", "UVM")

run_all_cox <- function(diseases){
  cox_ratios <- sapply(X = diseases, FUN = function(s) check_cox(paste0('../TCGA/', s,'.rData')))
  return (cox_ratios)
}

load("ratios.rData")
cox_ratios <- run_all_cox(diseases)
plot(cox_ratios, ratios[,1], ylab = 'Separation by median', xlab='Cox Model', main='Ratio of significant single genes')

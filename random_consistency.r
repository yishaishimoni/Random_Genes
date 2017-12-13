# Licensed Materials - Property of IBM 5748-XX8
# (C) Copyright IBM Corp. 1992, 1993 All Rights Reserved
# US Government Users Restricted Rights - Use, duplication or disclosure restricted by GSA ADP Schedule Contract with IBM Corp."

require(progress)
source('use_random_sets.r')
source('cleanData.r')

diseases <- c("ACC", "BLCA", "BRCA", "GBMLGG", "HNSC", "KIPAN", "KIRC",
              "KIRP", "LGG", "LIHC", "LUAD", "LUSC", "MESO",  
              "PAAD", "THYM", "UCEC", "UVM")

for(d in diseases){
  # check that the random signatures that are significant are consistent between training and testing
  load(paste0("../TCGA/",d,".rData"))
  dataObj <- cleanData(dataObj)
  print(paste("Number of samples in", d,":", nrow(dataObj$merged.dat)))
  
  M <- 50  # the number of splits into two halves
  N <- 100  # the number of random sets to take for each split
  pb <- progress_bar$new(total=M * N, width=80, format='[:bar] :percent :eta remaining')
  global_ps <- sapply(seq(M), FUN = function(y){
    # split the data randomly into two sets
    train_ind <- runif(nrow(dataObj$merged.dat)) > 0.5
    trainObj <- lapply(dataObj[3:4], FUN = function(x) x[train_ind,])
    testObj <- lapply(dataObj[3:4], FUN = function(x) x[!train_ind,])
  
    # repeat sampling random gene sets and testing their p-values on both data sets
    ps <- sapply(seq(N), FUN = function(z){
      pb$tick()
      robj <- singleRandomSet(trainObj$merged.dat, trainObj$Surv, returnGenes = TRUE, N = 64)
      p1 <- robj$pval
      genes <- robj$genes
      p2 <- singleRandomSet(testObj$merged.dat, testObj$Surv, genes=genes)
      return(c(p1, p2))
    })
    return(ps)
  })
  ps <- matrix(c(global_ps[1:N, ], global_ps[(N+1):(2*N), ]), N*M, 2)
  # plot the distribution
  plot(ps, log="xy")
  abline(h=0.05, v=0.05)
  # show a table of the proportion of significant sets in each half and in both
  contingency <- table(ps[, 1] < 0.05, ps[, 2] < 0.05)
  print(prop.table(contingency))
  ft <- fisher.test(contingency)
  print(ft)
  print(paste("p-value:", ft$p.value))
}

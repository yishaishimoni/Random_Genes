# Licensed Materials - Property of IBM 5748-XX8
# (C) Copyright IBM Corp. 1992, 1993 All Rights Reserved
# US Government Users Restricted Rights - Use, duplication or disclosure restricted by GSA ADP Schedule Contract with IBM Corp."

source('use_random_sets.r')
source('cleanData.r')
source('remove_meta_PCNA.r')

run_all <- function(N=50, reps=5000, remove_PCNA=F){
  message(paste0('Using N=', N, ' with remove_pcna=', remove_PCNA))
  message('-------------------------------')
  tcgaDir <- '../TCGA/'
  diseases <- c("ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "COADREAD", 
                "DLBC", "ESCA", "GBM", "GBMLGG", "HNSC", "KICH", "KIPAN", "KIRC",
                "KIRP", "LAML", "LGG", "LIHC", "LUAD", "LUSC", "MESO", "OV", 
                "PAAD", "PCPG", "PRAD", "READ", "SKCM", "TGCT", "THCA", 
                "THYM", "UCEC", "UCS", "UVM")
  filename <- paste0('randomGenesSummary_',N)
  if (remove_PCNA)
    filename <- paste0(filename, '_noPCNA')
  filename <- paste0(filename,'.pdf')
  pdf(filename); on.exit(dev.off())
  ratios <- vector(length=length(diseases))
  names(ratios) <- diseases
  for (d in diseases){
    # message('----------------------------')
    # message(paste("loading data for",d))
    load(paste0(tcgaDir,d,'.rData'))
    dataObj <- cleanData(dataObj)
    if (remove_PCNA)
      dataObj$merged.dat <- remove_meta_PCNA(dataObj)
    # message("running analysis")
    pvals <- multipleRandomSets(dataObj, d, reps, N)
    nulls <- pvals$null
    pvals <- pvals$pvals
    cutoff <- quantile(nulls, 0.05)
    ratio <- sum(pvals<cutoff)/length(pvals)
    ratios[d] <- ratio
    # significance <- get_pval(pvals, cutoff, ratio>0.05)
    # message(paste0(round(ratio*100),"% of pvals < 0.05"))
    ks_pval <- ks.test(nulls, pvals, alternative = 'two')$p.value
    chi_pval <- 1-pchisq(((ratio - 0.05) ^2) /0.05/0.95, df=1)
    Z = abs(ratio - 0.05) / sqrt(0.05*0.95/1000 + ratio*(1-ratio)/1000)
    prop_pval = 2*pnorm(-Z)
    # message(paste0('significance compared to null: ', ks_pval))
    message(paste(d, round(ratio*100), signif(prop_pval, digits=2), sep=' & '))
  }
  return (ratios)
}

get_pval <- function(pvals, cutoff, less){
  # given a cutoff, what is the probability that the proportion is actually 5%. 
  # we check this by bootsrapping
  ratios <- sapply(X=1:1000, FUN = function(x){
    b_pvals <- sample(x=pvals, size=length(pvals), replace=TRUE)
    return(sum(b_pvals<cutoff)/length(b_pvals))
  })
  # return(ratios)
  if(less){
    return((sum(ratios<0.05)+1) / length(ratios))
  } else {
    return((sum(ratios>0.05)+1) / length(ratios))
  }
}

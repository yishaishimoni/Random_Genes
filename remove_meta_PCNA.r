# Licensed Materials - Property of IBM 5748-XX8
# (C) Copyright IBM Corp. 1992, 1993 All Rights Reserved
# US Government Users Restricted Rights - Use, duplication or disclosure restricted by GSA ADP Schedule Contract with IBM Corp.

# remove meta_PCNA signature, in accordane with the original paper (DOI: 10.1371/journal.pcbi.1002240)
remove_meta_PCNA <- function(dataObj){
  meta_PCNA <- unlist(read.table('meta_PCNA.txt', stringsAsFactors=FALSE))
  meta_PCNA <- intersect(meta_PCNA, colnames(dataObj$merged.dat))
  PCNA_sig <- dataObj$merged.dat[ , meta_PCNA]
  # message(paste('removing PCNA signature using',ncol(PCNA_sig),'genes'))
  # for each sample, get the median of the meta_PCNA signature
  PCNA_sig <- apply(X = PCNA_sig, MARGIN = 1, FUN = function(x) median(x, na.rm = T))
  # for each gene, remove the linear fit of the PCNA signature across all samples
  dat <- apply(X = dataObj$merged.dat, MARGIN = 2, FUN = function(x, PCNA_sig){
    lm(x~PCNA_sig)$residuals + mean(x)
  }, PCNA_sig)
  return(dat)
}
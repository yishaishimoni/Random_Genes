# Licensed Materials - Property of IBM 5748-XX8
# (C) Copyright IBM Corp. 1992, 1993 All Rights Reserved
# US Government Users Restricted Rights - Use, duplication or disclosure restricted by GSA ADP Schedule Contract with IBM Corp.

#  automatically cluster each cancer type and check random genes within each cluster
require(phenoGraph)
source('use_random_sets.r')
source('cleanData.r')
source('remove_meta_PCNA.r')

random_subclasses <- TRUE
pdf_name = 'subclasses.pdf'
if(random_subclasses){
  pdf_name = 'random_subclasses.pdf'
}
runall <- function(){
  pdf(pdf_name); on.exit(dev.off())
  diseases <- c("ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "COADREAD", 
                "DLBC", "ESCA", "GBM", "GBMLGG", "HNSC", "KICH", "KIPAN", "KIRC",
                "KIRP", "LAML", "LGG", "LIHC", "LUAD", "LUSC", "MESO", "OV", 
                "PAAD", "PCPG", "PRAD", "READ", "SKCM", "TGCT", "THCA", 
                "THYM", "UCEC", "UCS", "UVM")
  
  tcgaDir <- '../TCGA/'
  ratios <- list()
  subclasses <- list()
  if(random_subclasses)
    load('subclasses.Rdata')
  for (disease in diseases){
    message(paste('loading', disease))
    load(paste0(tcgaDir,disease,'.rData'))
    dataObj <- cleanData(dataObj)
    if(random_subclasses){
      # define random assignments of the same sizes as the ones obtained by phenoclust
      C <- subclasses[[disease]]
      C$C <- sample(x=C$C, size=length(C$C), replace=FALSE)
    } else { 
      C <- phenoClust(t(dataObj$merged.dat), repeats = 10, k=10)
    }
    subclasses[[disease]] <- C
    # pvals <- list()
    ratios <- random_set_within_cluster(disease, ratios, C, dataObj)
    if(random_subclasses){
      save(ratios, file='random_subclassRatios.Rdata')
    } else {
      save(ratios, file='subclassRatios.Rdata')
      save(subclasses, file='subclasses.Rdata')
    }
    message('**************\n\n')
  }
}

runall()

# Licensed Materials - Property of IBM 5748-XX8
# (C) Copyright IBM Corp. 1992, 1993 All Rights Reserved
# US Government Users Restricted Rights - Use, duplication or disclosure restricted by GSA ADP Schedule Contract with IBM Corp.

source('use_random_sets.r')
source('cleanData.r')

# use the clinical stage information to create clusters and check their random bias

diseases <- c("ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "COADREAD", 
              "DLBC", "ESCA", "GBM", "GBMLGG", "HNSC", "KICH", "KIPAN", "KIRC",
              "KIRP", "LAML", "LGG", "LIHC", "LUAD", "LUSC", "MESO", "OV", 
              "PAAD", "PCPG", "PRAD", "READ", "SKCM", "TGCT", "THCA", 
              "THYM", "UCEC", "UCS", "UVM")

tcgaDir <- '../TCGA/'
ratios <- list()

pdf_name = 'clincal_subclasses.pdf'
runall <- function(){
  pdf(pdf_name); on.exit(dev.off())
  ratios <- list()
  subclasses <- list()
  for(disease in diseases){
    message(paste('loading', disease))
    load(paste0(tcgaDir,disease,'.rData'))
    dataObj <- cleanData(dataObj)
    grade_id <- grep(pattern = 'grade|stage', x = colnames(dataObj$clinical))
    if(length(grade_id) == 0) 
      next
    num_groups <- sapply(X=grade_id, FUN = function(i) length(levels(as.factor(dataObj$clinical[, i]))))
    keep <- num_groups > 2
    num_groups <- num_groups[keep]
    grade_id <- grade_id[keep]
    if(length(grade_id) == 0) 
      next
    grade_id <- grade_id[which.min(num_groups)]
    C <- list()
    C$C <- as.numeric(addNA(as.factor(dataObj$clinical[rownames(dataObj$merged.dat), grade_id])))
    ratios <- random_set_within_cluster(disease, ratios, C, dataObj)
    save(ratios, file='clinical_subclassRatios.Rdata')
    save(subclasses, file='clinical_subclasses.Rdata')
    message('**************\n\n')
  }
}

runall()
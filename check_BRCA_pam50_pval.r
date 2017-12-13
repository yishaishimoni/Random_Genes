# Licensed Materials - Property of IBM 5748-XX8
# (C) Copyright IBM Corp. 1992, 1993 All Rights Reserved
# US Government Users Restricted Rights - Use, duplication or disclosure restricted by GSA ADP Schedule Contract with IBM Corp.

source('use_random_sets.r')
source('cleanData.r')
source('remove_meta_PCNA.r')

tcgaDir <- '../TCGA/'

load(paste0(tcgaDir,'BRCA.rData'))
dataObj <- cleanData(dataObj)

pam50_centroids <- read.delim("D:/random_Genes/pam50_centroids.txt", row.names=1, stringsAsFactors=FALSE)
pam50_genes <- rownames(pam50_centroids)

I <- intersect(rownames(pam50_centroids), colnames(dataObj$merged.dat))
pam50_brca <- dataObj$merged.dat[, I]
pam50_centroids <- pam50_centroids[I, ]
singleRandomSet(pam50_brca, dataObj$Surv, toPlot = TRUE, N=48, replace=FALSE)

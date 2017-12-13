# Licensed Materials - Property of IBM 5748-XX8
# (C) Copyright IBM Corp. 1992, 1993 All Rights Reserved
# US Government Users Restricted Rights - Use, duplication or disclosure restricted by GSA ADP Schedule Contract with IBM Corp.

# identify the BRCA subclasses by pam50 and compare them with the phenoClust results
require(survival)
load('../TCGA/BRCA.rData')
source('cleanData.r')
dataObj <- cleanData(dataObj)
pam50_centroids <- read.delim("D:/random_Genes/pam50_centroids.txt", row.names=1, stringsAsFactors=FALSE)

I <- intersect(rownames(pam50_centroids), colnames(dataObj$merged.dat))
pam50_brca <- dataObj$merged.dat[, I]
pam50_centroids <- pam50_centroids[I, ]

require(amap)
D <- cor(x=t(pam50_brca), y=pam50_centroids)
C <- apply(X = D, MARGIN = 1, FUN = which.max)

source('use_random_sets.r')
remove_PCNA <- F
ratios <- list()
for (cluster in 1:max(C)){
  message(paste('cheking cluster', cluster))
  I <- C==cluster
  dataObj_C <- dataObj
  dataObj_C$clinical <- dataObj_C$clinical[I]
  dataObj_C$merged.dat <- dataObj_C$merged.dat[I, ]
  dataObj_C$Surv <- dataObj_C$Surv[I, ]
  if (remove_PCNA){
    dataObj_C$merged.data <- remove_meta_PCNA(dataObj_C)
  }
  pvals <- multipleRandomSets(dataObj=dataObj_C, 
                              d=paste0('BRCA (', colnames(pam50_centroids)[cluster],')'))
  ratios[[colnames(pam50_centroids)[cluster]]] <- sum(pvals<0.05)/length(pvals)
}
save(D, C, ratios, file = 'pam50ratios.Rdata')

pdf('pam50_vs_phenoclust.pdf')
barplot(unlist(ratios))

# ratios are good - they are all close to 0.05, indicating that there are no obvious
# additional subclasses
load('subclasses.Rdata')
C_table <- table(subclasses$BRCA$C, C)
colnames(C_table) <- names(pam50_centroids)
print(C_table)
barplot(C_table, col=1:max(subclasses$BRCA$C))
legend(x="topleft", fill=1:max(subclasses$BRCA$C), legend = paste('class', 1:max(subclasses$BRCA$C)))
# It can be seen that basal were identified in class 5, HER2 by class 2, luminal B mostly in class 1, and 
# Luminal A mostly in classes 3 and 4. 

barplot(t(C_table), col=1:max(C))
legend(x="topright", fill=1:max(C), legend=colnames(pam50_centroids), cex=0.9)
# conversely, cluster 1 is dominated by Luminal B, class 2 by HER2, classes 3 and 4 by Luimnal A, 
# class 5 by Basal, while class 6 is almost evenly split between Luminal A and B.
dev.off()

# Licensed Materials - Property of IBM 5748-XX8
# (C) Copyright IBM Corp. 1992, 1993 All Rights Reserved
# US Government Users Restricted Rights - Use, duplication or disclosure restricted by GSA ADP Schedule Contract with IBM Corp.

# plot ratios
require(gplots)
require(extrafont)
plotRatios <- function(filename){
  load(filename)
  curPar <- par(no.readonly = T); on.exit(par(curPar))
  mar = par("mar") + c(-0.8, -0.2, 3, 24)
  par(mar=mar)
  matplot(x=2^(1:ncol(ratios)-1), y=t(ratios), type="l", log="x", col=c(1:6,8), lty=1:5, lwd=2,
          xlab="No. of random genes", ylab="Proportion of significant sets", ylim=c(0,1))
  # rect(1,0,2^(ncol(ratios)-1),0.05,col = '#55555555', border=F )
  
  par(new=T, mar=c(0,0,0,0))
  plot.new()
  legend(x="top",ncol=7, col=c(1:6,8), lty=1:5, legend=rownames(ratios), bty = "n", lwd=2, cex=0.9)

  heatmap.2(x = ratios, Colv = NA, scale = "none", trace="n", dendrogram = "none", labCol = 2^(1:ncol(ratios)-1),
            col=colorpanel(n=16,low = "white",high = "cornflowerblue"))
}

plotBothRatios <- function(){
  pdf('signifRatiosPerTumor.pdf', width=10, height = 5); on.exit(dev.off())
  plotRatios('ratios.rData')
  plotRatios('ratios_noPCNA.Rdata')
}

plotBothRatios()

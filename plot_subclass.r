# Licensed Materials - Property of IBM 5748-XX8
# (C) Copyright IBM Corp. 1992, 1993 All Rights Reserved
# US Government Users Restricted Rights - Use, duplication or disclosure restricted by GSA ADP Schedule Contract with IBM Corp.

require(extrafont)
# plot the results of the subclass analysis
plot_subclass <- function (prefix='', main='PhenoClust') {
  # prefix <- 'random_'; main <- 'random'
  # prefix <- 'clinical_'; main <- 'Grade-based'
  width=3.25 * 12 / 8
  height=4.5 * 12 / 8
  load(paste0(prefix,'subclassRatios.Rdata'))
  # load('random_subclassRatios.Rdata')
  pdf(paste0(prefix, 'subclass_improve.pdf'), width=width, height=height, family = "Arial"); on.exit(dev.off())
  par(mar=par("mar") + c(0, 1.6, 0, -0.6))
  plot(range(pmax(unlist(ratios), 0.01)),range(pmax(unlist(ratios), 0.01)),log="xy", type="l",
       main='Ratio of significant sets', xlab = 'without correction', ylab='after correction by PCNA')
  for (i in seq(1,length(ratios),2)){
    points(unlist(ratios[[i]]),unlist(ratios[[i+1]]), xlab="ratio",ylab="PCNA corrected")
  }
  abline(h=0.05, v=0.05)
  
  subRatios <- ratios
  load('ratios.rData')
  for(type in rownames(ratios)){
    if(!(type %in% names(subRatios))){
      subRatios[[type]] <- list(c(NaN))
      subRatios[[paste0(type,'noPCNA')]] <- list(c(NaN))
    }
  }

  rs <- NULL
  group <- NULL
  for (group_name in rownames(ratios)){
    temp <- unlist(subRatios[group_name])
    if(all(is.na(temp))){
      rs <- c(rs, 2)
      group <- c(group, group_name)
    } else {
      rs <- c(rs, temp)
      group <- c(group, rep(group_name, length(temp)))
    }
  }
  stripchart(rep(1, length(group)) ~ group, vertical = FALSE, pch=16, col="white", 
             las=1, xlab='Proportion of significant random sets', xlim=c(0, 1), 
             main=paste(main, 'clusters'))
  stripchart(rs ~ group, vertical = FALSE, pch=16, col=1:nrow(ratios), add=TRUE)
  for (i in 1:nrow(ratios)){
    abline(h=i, col='lightgrey', lty=3)
    lines(y=c(i-0.4, i+0.4), x=c(ratios[i, 7], ratios[i, 7]), col=8, lwd=2)
    # lines(y=c(i-0.4, i+0.4), x=c(ratios[i, 7], ratios[i, 7]), col=i, lwd=2)
  }
  abline(v=0.05, col='grey', lwd=1, lty=2)
  stripchart(rs ~ group, vertical = FALSE, pch=16, col=1:nrow(ratios), add=TRUE)
  stripchart(rs ~ group, vertical = FALSE, pch=1, add=TRUE)
  legend(x = 0.5, y=31, legend=c('original dataset','sub-cluster'), pch = c(NA, 16), lwd=c(2, NA), col=c(8, 1), )
  
  plot(density(na.omit(unlist(subRatios[1:nrow(ratios)*2 - 1]))))
  
  b <- barplot(ratios[ , 7], las=2, col=0, border = NA, horiz = TRUE, xlim=c(0,1))
  d <- (b[2]-b[1])/2 * 0.6
  # usr <- par()$usr
  # usr[1:2] <- c(0,nrow(ratios))
  # par(usr=usr)
  abline(v=0.05, col='grey', lwd=2)
  for (i in 1:nrow(ratios)){
    group_name <- rownames(ratios)[i]
    r <- ratios[i, 7]
    rs <- unlist(subRatios[group_name])
    xs <- seq(b[i]-d, b[i]+d, length.out = length(rs)+1)
    arrows_xs = seq(b[i]+d, b[i]-d, length.out = length(rs))
    if(length(rs)>0){
      arrows(y0 = arrows_xs, x0=r, x1=rs, col=1, length = 0.05, lwd=2)
      arrows(y0 = arrows_xs, x0=r, x1=rs, col=((i-1) %% 7) + 2, length = 0.05, lwd=1)
      # rect(xleft=xs[1:length(rs)], ybottom=rs, xright=xs[-c(1)], ytop=r, col=((i-1) %% 7) + 2)
    }
    lines(y=c(min(xs), max(xs)), x=c(r, r), col='black', lwd=3)
  }
}

plot_subclass()
plot_subclass(prefix='random_', main='random')
plot_subclass(prefix='clinical_', main='Grade-based')


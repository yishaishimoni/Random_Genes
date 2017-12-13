# Licensed Materials - Property of IBM 5748-XX8
# (C) Copyright IBM Corp. 1992, 1993 All Rights Reserved
# US Government Users Restricted Rights - Use, duplication or disclosure restricted by GSA ADP Schedule Contract with IBM Corp.

# check if random genes predict survival

source('run_all.r')

run_once <- function(remove_PCNA=TRUE){
	ratios <- sapply(X = 2^c(0:10), FUN = function(N) run_all(N=N, reps = 5000, remove_PCNA = remove_PCNA))
	filename <- 'ratios'
	if (remove_PCNA)
		filename <- paste0(filename,'_noPCNA')
	filename <- paste0(filename,'.Rdata') 
	save(ratios, file=filename)
}

run_once()
run_once(FALSE)
# Licensed Materials - Property of IBM 5748-XX8
# (C) Copyright IBM Corp. 1992, 1993 All Rights Reserved
# US Government Users Restricted Rights - Use, duplication or disclosure restricted by GSA ADP Schedule Contract with IBM Corp.

# fix cases where the merge was broken
if (nrow(dataObj$merged.dat)==0){
  shortPatientIDs <- substr(colnames(dataObj$dat),start = 1, stop = 12)
  patientsWithClinData <- intersect(rownames(dataObj$clinical), shortPatientIDs)
  IdataObj <- match(patientsWithClinData, shortPatientIDs)
  bcr <- patientsWithClinData
  Iclinical <- match(patientsWithClinData, bcr)
  status <- as.numeric(dataObj$clinical[ ,"vitalstatus"])
  OS <- as.numeric(dataObj$clinical[ ,"daystodeath"])
  OS[status==0] <- as.numeric(dataObj$clinical[status==0,"daystolastfollowup"])
  OS[is.na(OS)] <- as.numeric(dataObj$clinical[is.na(OS),"daystolastknownalive"])
  status <- status[Iclinical]
  OS <- OS[Iclinical]
  dataObj$merged.dat <- as.data.frame(cbind(status,OS,t(dataObj$dat[ ,IdataObj])))
  dataObj$merged.dat <- cbind(bcr, dataObj$merged.dat)
}

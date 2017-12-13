# Licensed Materials - Property of IBM 5748-XX8
# (C) Copyright IBM Corp. 1992, 1993 All Rights Reserved
# US Government Users Restricted Rights - Use, duplication or disclosure restricted by GSA ADP Schedule Contract with IBM Corp.

# run this on a linux machine
# load TCGA data
require(TCGA2STAT)
source('getTCGA.r') # I had to use a local updated version 
Sys.setenv(TAR="c:/Cygwin64/bin/tar", R_GZIPCMD="c:/Cygwin64/bin/gzip")
tcgaDir <- '../TCGA/'

diseases <- c("ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "COADREAD", 
              "DLBC", "ESCA", "GBM", "GBMLGG", "HNSC", "KICH", "KIPAN", "KIRC",
              "KIRP", "LAML", "LGG", "LIHC", "LUAD", "LUSC", "MESO", "OV", 
              "PAAD", "PCPG", "PRAD", "READ", "SKCM", "TGCT", "THCA", 
              "THYM", "UCEC", "UCS", "UVM")

for (d in diseases){
  message(paste("importing",d))
  dataObj <- getTCGA(disease = d, clinical = T)
  save(dataObj, file=paste0(tcgaDir,d,".rData"))
}

# Abbreviation 	Cancer Type
# ACC	Adrenocortical carcinoma
# BLCA	Bladder Urothelial Carcinoma
# BRCA	Breast invasive carcinoma
# CESC	Cervical squamous cell carcinoma and endocervical adenocarcinoma
# CHOL	Cholangiocarcinoma
# CNTL	Controls
# COAD	Colon adenocarcinoma
# DLBC	Lymphoid Neoplasm Diffuse Large B-cell Lymphoma
# ESCA	Esophageal carcinoma
# GBM	Glioblastoma multiforme
# HNSC	Head and Neck squamous cell carcinoma
# KICH	Kidney Chromophobe
# KIRC	Kidney renal clear cell carcinoma
# KIRP	Kidney renal papillary cell carcinoma
# LAML	Acute Myeloid Leukemia
# LCML	Chronic Myelogenous Leukemia
# LGG	Brain Lower Grade Glioma
# LIHC	Liver hepatocellular carcinoma
# LUAD	Lung adenocarcinoma
# LUSC	Lung squamous cell carcinoma
# MESO	Mesothelioma
# MISC	Miscellaneous
# OV	Ovarian serous cystadenocarcinoma
# PAAD	Pancreatic adenocarcinoma
# PCPG	Pheochromocytoma and Paraganglioma
# PRAD	Prostate adenocarcinoma
# READ	Rectum adenocarcinoma
# SARC	Sarcoma
# SKCM	Skin Cutaneous Melanoma
# STAD	Stomach adenocarcinoma
# TGCT	Testicular Germ Cell Tumors
# THCA	Thyroid carcinoma
# THYM	Thymoma
# UCEC	Uterine Corpus Endometrial Carcinoma
# UCS	Uterine Carcinosarcoma
# UVM	Uveal Melanoma
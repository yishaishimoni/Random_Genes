# Code and data related to manuscript

This is an updated version of the code that was used for the analysis of TCGA gene-expression (RNA-seq) data, 
checking what proportion of random gene sets are significantly predictive of cancer survival.
See the full paper 
[Association between expression of random gene sets and survival is evident in multiple cancer types and may be explained by sub-classification](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006026)
in PLoS Computational Biology. 

## Reproduce figures
The `.Rdata` files are the saved results of the original analysis.
If you only wish to recreate and perhaps manupulate the figures, 
all you need is to run the scripts `plotRatios.r` and `plot_subclass.r'.
This will create multiple pdf files with the results. 

## Rerun analysis
It is possible to run the whole analysis from scratch using this code.
However the TCGA data is a big dataset (and so cannot fit on the git repo) 
and it is not mine to share (and so should not sit on the git repo).
You can download the data using the script 'load_data.r'.

I will point out that the original analysis was done using a deprecated package and a deprecated TCGA interface.
Therefore, the exact same data cannot be replicated. 
The `load_data.r` script downloads the current data and allows runnign the same analysis on this new data.
Downloading the data take quite some time, and includes both downloading and preprocessing steps.
This script will populate with new .Rdata files in the TCGA folder that include the expression and survival data of all cancer types.

To re-create the analysis do the following:
- Run the `testRandomGenes.r` script. This will create the ratio files with and without PCNA correction, and for different sizes of random sets.
  This can take quite a long time to run.
- Run the script `plotRatios.r` to create plots of the ratios that are saved in the ratios files.
- Run `random_consistency.r` to check that the random bias is consistent across the data.
- Run the script `use_subclasses.r` that automatically identifies subclasses, checks the significance ratios within each subclass, 
  then checks the significance ratios for random subclasses of the same size. 
  It also creates multiple pdf files with plots of the resulting distributions of the p-values, and saves the subclasses and their results to files.
  This can take a long time to run.
- Run `clinical_variables_as_subclasses.r` to check the result of sub-classification by grade.
- To create the subclass plots run `plot_subclass.r`.

### Requirements
The code requires the packages `survival`, `amap`, `progress`, `gplots`, and `extrafonts`.
You will also need to install `remotes`, and `BiocManager` to install additional dependencies using 
```R
remotes::install_github('yishaishimoni/phenoClust')
BiocManager::install('TCGAbiolinks')
```

### Install using Dockerfile
Alternatively, I created a [Dockerfile](https://github.com/yishaishimoni/dockers/tree/main/random_genes/Dockerfile) that will build an environment with Rstudio, all the packages, and the code.
You will need to build and start the docker container as explained in the [dockerfile's repository](https://github.com/yishaishimoni/dockers/tree/main/random_genes).

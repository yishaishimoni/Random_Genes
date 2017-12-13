## Code and data related to manuscript

This is the code that was used for the analysis of TCGA gene-expression (RNA-seq) data, 
checking what proportion of random gene sets are significantly predictive of cancer survival.

The `.Rdata` files are the saved results of this analysis.
However, it is possible to run the whole analysis from scratch using this code.
The main script is `testRandomGenes.r`, which creates the significance ratios for all cancer types (with or without PCNA correction), 
and also plots the resulting p-value distributions to pdf files.

Before running the analysis, make sure to have the TCGA data downloaded into `../TCGA`.
This can be done using `load_data.r`. Please note that I had to tweak the code a bit for it to run on a windows machine 
(specifically, I had to change the pointers to `gzip` and `tar`)
This is a big dataset (and so cannot fit on the git repo) 
and it is not mine to share (and so cannot sit on the git repo).

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
The scripts require the package `PhenoGraph`, which can be installed directly from git by installing the `devtools` package and running
```R
> devtools::install_github('yishaishimoni/phenoClust')
```

The code also requires the packages `survival`, `TCGA2STAT`, `amap`, `progress`, `gplots`, and `extrafonts`.

To properly download the data requires both `gzip` and `tar` (at least in the version I had for `TCGA2STAT`).
These have to be installed separately (e.g. using cygwin) and pointed to in the `load_data.r` script.
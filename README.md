# asmbPLS: Predicting and Classfying Patient Phenotypes with Multi-omics Data

Runzhi Zhang, Susmita Datta

## Description
Adaptive Sparse Multi-block Partial Least Square, a supervised algorithm, is an extension of the Sparse Multi-block Partial Least Square, which allows different quantiles to be used in different blocks of different partial least square components to decide the proportion of features to be retained. The best combinations of quantiles can be chosen from a set of user-defined quantiles combinations by cross-validation. By doing this, it enables us to do the feature selection for different blocks, and the selected features can then be further used to predict the outcome. For example, in biomedical applications, clinical covariates plus different types of omics data such as microbiome, metabolome, mRNA data, methylation data, copy number variation data might be predictive for patients outcome such as survival time or response to therapy. Different types of data could be put in different blocks and along with survival time to fit the model. The fitted model can then be used to predict the survival for the new samples with the corresponding clinical covariates and omics data. 

In addition, Adaptive Sparse Multi-block Partial Least Square Discriminant Analysis is also included, which extends Adaptive Sparse Multi-block Partial Least Square for classifying the categorical outcome.

## R package: asmbPLS (Version: 1.0.1)

## R package installation
```
install_github("RunzhiZ/asmbPLS")
```
or
```
install.packages("asmbPLS")
```
R cran link: https://cran.r-project.org/web/packages/asmbPLS/index.html

### Common errors for MAC users:
* **Error 1**:
```
ld: library not found for -lgfortran
```
Solution for **error 1**: install the required tools https://mac.r-project.org/tools/


* **Error 2**:
```
clang: error: unsupported option '-fopenmp'
```
Possible solution for **error 2**: https://stackoverflow.com/questions/43555410/enable-openmp-support-in-clang-in-mac-os-x-sierra-mojave

## Installation error report
If you have more errors installing the R package, please report to runzhi.zhang@ufl.edu

## Tutorial
[Click here to view the tutorial for the R package](https://rpubs.com/spencer886/1026803)

## References
* [Zhang R, Datta S: asmbPLS: Adaptive Sparse Multi-block Partial Least Square for Survival Prediction using Multi-Omics Data. bioRxiv 2023:2023.2004.2003.535442.](https://www.biorxiv.org/content/10.1101/2023.04.03.535442v1)
* [Zhang R, Datta S: Adaptive Sparse Multi-Block PLS Discriminant Analysis: An Integrative Method for Identifying Key Biomarkers from Multi-Omics Data. Genes 2023, 14(5), 961.](https://www.mdpi.com/2073-4425/14/5/961/htm)

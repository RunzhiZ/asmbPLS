# asmbPLS: Predicting and Classfying Patient Phenotypes with Multi-omics Data

Runzhi Zhang, Susmita Datta

## Description
Adaptive Sparse Multi-block Partial Least Square (asmbPLS), a supervised algorithm, is an extension of the smbPLS, which allows different quantiles to be used in different blocks of different PLS components to decide the proportion of features to be retained. The best combinations of quantiles can be chosen from a set of user-defined quantiles combinations by cross-validation. By doing this, asmbPLS enables us to do the feature selection for differnt blocks, and the selected features can then be further used to predict the outcome. For example, in biomedical applications, clinical covariates plus different types of omics data such as microbiome, metabolome, mRNA data, methylation data, CNV data might be predictive for patients outcome such as survival timeor response to therapy. Different types of data could be put in differnt blocks and along with survival time to fit the asmbPLS model. The fitted model can then be used to predict the survival for the new samples with the corresponding clinical covariates and omics data.

In addition, Adaptive Sparse Multi-block Partial Least Square Discriminant Analysis (asmbPLS-DA) is also included, which extends asmbPLS for classifying the categorical outcome.

## R package installation
```
devtools::install_github("RunzhiZ/asmbPLS")
```
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

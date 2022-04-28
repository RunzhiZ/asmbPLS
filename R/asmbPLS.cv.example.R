#' Example data for asmbPLS CV
#'
#' Simulated data for asmbPLS CV. 
#'
#' @docType data
#'
#' @usage data(asmbPLS.cv.example)
#'
#' @format A list including 6 components:  
#' 1) X.matrix, a matrix with 100 samples (rows) and 400 features (columns, 1-200 are microbial taxa, 201-400 are metabolites); 
#' 2) Y.matrix, a matrix with 100 samples (rows) and 1 column (log-transformed survival time); 
#' 3) Y.indicator, a vector containing the event indicator for each sample
#' 4) X.dim, dimension of the two blocks in X.matrix; 
#' 5) quantile.comn.table, user-defined quantile combinations used for CV; 
#' 6) PLS.comp, selected number of PLS components
#' 
#' @keywords datasets
#'
"asmbPLS.cv.example"
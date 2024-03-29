#' Example data for asmbPLS-DA algorithm
#'
#' Simulated data for asmbPLS-DA. 
#'
#' @docType data
#'
#' @usage data(asmbPLSDA.example)
#'
#' @format A list including 8 components:  
#' 
#' 1) \code{X.matrix}, a matrix with 100 samples (rows) and 400 features, features 1-200 are from block 1 and features 201-400 are from block 2; 
#' 
#' 2) \code{X.matrix.new}, a matrix to be predicted with 100 samples (rows) and 400 features, features 1-200 are from block 1 and features 201-400 are from block 2;  
#' 
#' 3) \code{Y.matrix.binary}, a matrix with 100 samples (rows) and 1 column; 
#' 
#' 4) \code{Y.matrix.morethan2levels}, a matrix with 100 samples (rows) and 3 columns (3 levels);
#' 
#' 5) \code{X.dim}, dimension of the two blocks in X.matrix; 
#' 
#' 6) \code{PLS.comp}, selected number of PLS components; 
#' 
#' 7) \code{quantile.comb}, selected quantile combinations;
#' 
#' 8) \code{quantile.comb.table.cv}, pre-defined quantile combinations for cross validation.
#' 
#' @keywords datasets
#'
"asmbPLSDA.example"
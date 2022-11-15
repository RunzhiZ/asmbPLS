#' Example data for asmbPLS algorithm
#'
#' Simulated data for asmbPLS. 
#'
#' @docType data
#'
#' @usage data(asmbPLS.example)
#'
#' @format A list including 8 components:  
#' 
#' 1) \code{X.matrix}, a matrix with 100 samples (rows) and 400 features (columns, 1-200 are microbial taxa, 201-400 are metabolites); 
#' 
#' 2) \code{X.matrix.new}, a matrix to be predicted with 100 samples (rows) and 400 features (columns, 1-200 are microbial taxa, 201-400 are metabolites); 
#' 
#' 3) \code{Y.matrix}, a matrix with 100 samples (rows) and 1 column (log-transformed survival time); 
#' 
#' 4) \code{X.dim}, dimension of the two blocks in X.matrix; 
#' 
#' 5) \code{PLS.comp}, selected number of PLS components;
#' 
#' 6) \code{quantile.comb}, selected quantile combinations; 
#' 
#' 7) \code{quantile.comb.table.cv}, pre-defined quantile combinations for cross validaiton;
#' 
#' 8) \code{Y.indicator}, a vector containing the event indicator for each sample.
#' 
#' @keywords datasets
#'
"asmbPLS.example"
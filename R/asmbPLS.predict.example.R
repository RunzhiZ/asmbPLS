#' Example data for asmbPLS predict
#'
#' Simulated data for asmbPLS predict. 
#'
#' @docType data
#'
#' @usage data(asmbPLS.predict.example)
#'
#' @format A list including 6 components:  
#' 1) X.matrix, a matrix with 100 samples (rows) and 400 features (columns, 1-200 are microbial taxa, 201-400 are metabolites); 
#' 2) X.matrix.new a matrix to be predicted with 100 samples (rows) and 400 features (columns, 1-200 are microbial taxa, 201-400 are metabolites); 
#' 3) Y.matrix, a matrix with 100 samples (rows) and 1 column (log-transformed survival time); 
#' 4) X.dim, dimension of the two blocks in X.matrix; 
#' 5) quantile.comn, selected quantile combinations; 
#' 6) PLS.comp, selected number of PLS components
#' 
#' @keywords datasets
#'
"asmbPLS.predict.example"
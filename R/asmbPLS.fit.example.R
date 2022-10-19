#' Example data for asmbPLS fit
#'
#' Simulated data for asmbPLS fit. 
#'
#' @docType data
#'
#' @usage data(asmbPLS.fit.example)
#'
#' @format A list including 5 components:  
#' 1) X.matrix, a matrix with 100 samples (rows) and 400 features (columns, 1-200 are microbial taxa, 201-400 are metabolites); 
#' 2) Y.matrix, a matrix with 100 samples (rows) and 1 column (log-transformed survival time); 
#' 3) X.dim, dimension of the two blocks in X.matrix; 
#' 4) quantile.comp, selected quantile combinations; 
#' 5) PLS.comp, selected number of PLS components
#' 
#' @keywords datasets
#'
"asmbPLS.fit.example"
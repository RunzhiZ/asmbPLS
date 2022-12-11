#' asmbPLS-DA for block-structured data
#'
#' Function to fit the adaptive sparse multi-block partial least square 
#' discriminant analysis (asmbPLS-DA) model with several explanatory blocks 
#' (X_1, ..., X_B) as our predictors to explain the categorical outcome Y.
#' 
#' @param X.matrix Predictors matrix. Samples in rows, variables in columns
#' @param Y.matrix Outcome matrix. Samples in rows, this is a matrix with one 
#' column (binary) or multiple columns (more than 2 levels, dummy variables). 
#' @param PLS.comp Number of PLS components in asmbPLS-DA.
#' @param X.dim A vector containing the number of predictors in each block 
#' (ordered).
#' @param quantile.comb A matrix containing quantile combinations used for 
#' different PLS components, whose row number equals to the number of PLS 
#' components used, column number equals to the number of blocks.
#' @param outcome.type The type of the outcome Y. "\code{binary}" for binary 
#' outcome, and "\code{multiclass}" for categorical outcome with more than 2 
#' levels.
#' @param center A logical value indicating whether weighted mean center should 
#' be implemented for X.matrix and Y.matrix. The default is TRUE.
#' @param scale  A logical value indicating whether scale should be 
#' implemented for X.matrix. The default is TRUE.
#' @param maxiter A integer indicating the maximum number of iteration. The
#' default number is 100.
#' 
#' @return 
#' \code{asmbPLSDA.fit} returns a list containing the following components:
#' \item{X_dim}{A vector containing the number of predictors in each block.}
#' \item{X_weight}{A list containing the weights of predictors for different 
#' blocks in different PLS components.}
#' \item{X_score}{A list containing the scores of samples in different blocks
#' in different PLS components.}
#' \item{X_loading}{A list containing the loadings of predictors for different
#' blocks in different PLS components.}
#' \item{X_super_weight}{A matrix containing the super weights of different
#' blocks for different PLS components.}
#' \item{X_super_score}{A matrix containing the super scores of samples for
#' different PLS components.}
#' \item{Y_weight}{A matrix containing the weights of outcome for different 
#' PLS components.}
#' \item{Y_score}{A matrix containing the scores of outcome for different 
#' PLS components.}
#' \item{X_col_mean}{A matrix containing the weighted mean of each predictor 
#' for scaling.}
#' \item{Y_col_mean}{The weighted mean of outcome matrix for scaling.}
#' \item{X_col_sd}{A matrix containing the standard deviation (sd) of each 
#' predictor for scaling. sd for predictors with sd = 0 will be changed to 1.}
#' \item{center}{A logical value indicating whether weighted mean center is
#' implemented for X.matrix.}
#' \item{scale}{A logical value indicating whether scale is implemented for 
#' X.matrix.}
#' \item{Outcome_type}{The type of the outcome Y. "\code{binary}" for binary 
#' outcome, and "\code{multiclass}" for categorical outcome with more than 2
#' levels.}
#' \item{Y_group}{Original Y.matrix.}

#' @examples
#' ## Use the example dataset
#' data(asmbPLSDA.example)
#' X.matrix = asmbPLSDA.example$X.matrix
#' Y.matrix.binary = asmbPLSDA.example$Y.matrix.binary
#' Y.matrix.multiclass = asmbPLSDA.example$Y.matrix.morethan2levels
#' X.dim = asmbPLSDA.example$X.dim
#' PLS.comp = asmbPLSDA.example$PLS.comp
#' quantile.comb = asmbPLSDA.example$quantile.comb
#'  
#' ## asmbPLSDA fit for binary outcome
#' asmbPLSDA.fit.binary <- asmbPLSDA.fit(X.matrix = X.matrix, 
#'                                       Y.matrix = Y.matrix.binary, 
#'                                       PLS.comp = PLS.comp, 
#'                                       X.dim = X.dim, 
#'                                       quantile.comb = quantile.comb,
#'                                       outcome.type = "binary")
#' 
#' ## asmbPLSDA fit for categorical outcome with more than 2 levels
#' asmbPLSDA.fit.multiclass <- asmbPLSDA.fit(X.matrix = X.matrix, 
#'                                           Y.matrix = Y.matrix.multiclass, 
#'                                           PLS.comp = PLS.comp, 
#'                                           X.dim = X.dim, 
#'                                           quantile.comb = quantile.comb,
#'                                           outcome.type = "multiclass")
#' 
#' @export
#' @useDynLib asmbPLS, .registration=TRUE
#' @importFrom Rcpp sourceCpp

asmbPLSDA.fit <- function(X.matrix, Y.matrix, PLS.comp, X.dim, quantile.comb, outcome.type, center = TRUE, scale = TRUE, maxiter = 100){
  stopifnot(!missing(X.matrix),
            !missing(Y.matrix),
            !missing(PLS.comp),
            !missing(X.dim),
            !missing(quantile.comb),
            !missing(outcome.type),
            is.matrix(X.matrix), 
            is.matrix(Y.matrix),
            is.matrix(quantile.comb))
  return(asmbPLSDA_fit(X.matrix, Y.matrix, PLS.comp, X.dim, quantile.comb, outcome.type, center, scale, maxiter))
}
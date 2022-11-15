#' Cross-validation for asmbPLS to find the best combinations of quantiles for prediction
#'
#' Function to find the best combinations of quantiles used for prediction via
#' cross-validation. Usually should be conducted before 
#' \code{\link[asmbPLS]{asmbPLS.fit}} to obtain the quantile combinations.
#' 
#' @param X.matrix Predictors matrix. Samples in rows, variables in columns.
#' @param Y.matrix Outcome matrix. Samples in rows, this is a matrix with one 
#' column (continuous variable). The outcome could be imputed survival time or 
#' other types of continuous outcome. For survival time with right-censored 
#' survival time and event indicator, the right censored time could be imputed 
#' by \code{\link{meanimp}}.
#' @param PLS.comp Number of PLS components in asmbPLS.
#' @param X.dim A vector containing the number of predictors in each block 
#' (ordered).
#' @param quantile.comb.table A matrix containing user-defined quantile 
#' combinations used for CV, whose column number equals to the 
#' number of blocks.
#' @param Y.indicator A vector containing the event indicator for each sample, 
#' whose length is equal to the number of samples. This vector allows the ratio 
#' of observed/unobserved to be the same in the training set and validation set. 
#' Observed = 1, and unobserved = 0. If other types of outcome data rather than 
#' survival outcome is used, you can use a vector with all components = 1 
#' instead.
#' @param k The number of folds of CV procedure. The default is 5.
#' @param only.observe Whether only observed samples in the validation set 
#' should be used for calculating the MSE for CV. The default is TRUE.
#' @param seed An integer given by user to obtain reproducible results. The
#' default is 1.
#' 
#' @return 
#' \code{asmbPLS.cv} returns a list containing the following components:
#' \item{quantile_table_CV}{A matrix containing the selected quantile 
#' combination and the corresponding MSE of CV for each PLS component.}
#' \item{CV_results}{A list containing the details of the CV results for each PLS 
#' component.}
#' \item{CV_index}{A list containing the validation and training index for each 
#' cross validation fold.}
#' 
#' @examples
#' ## Use the example dataset
#' data(asmbPLS.example)
#' X.matrix = asmbPLS.example$X.matrix
#' Y.matrix = asmbPLS.example$Y.matrix
#' PLS.comp = asmbPLS.example$PLS.comp
#' X.dim = asmbPLS.example$X.dim
#' quantile.comb.table.cv = asmbPLS.example$quantile.comb.table.cv
#' Y.indicator = asmbPLS.example$Y.indicator
#' 
#' ## cv to find the best quantile combinations for model fitting
#' cv.results <- asmbPLS.cv(X.matrix = X.matrix, 
#'                          Y.matrix = Y.matrix, 
#'                          PLS.comp = PLS.comp, 
#'                          X.dim = X.dim, 
#'                          quantile.comb.table = quantile.comb.table.cv, 
#'                          Y.indicator = Y.indicator, 
#'                          k = 5, 
#'                          only.observe = TRUE, 
#'                          seed = 123)
#' quantile.comb <- cv.results$quantile_table_CV[,1:length(X.dim)]
#'  
#' ## asmbPLS fit
#' asmbPLS.results <- asmbPLS.fit(X.matrix = X.matrix, 
#'                                Y.matrix = Y.matrix, 
#'                                PLS.comp = PLS.comp, 
#'                                X.dim = X.dim, 
#'                                quantile.comb = quantile.comb)
#' 
#' @export
#' @useDynLib asmbPLS, .registration=TRUE
#' @importFrom Rcpp sourceCpp

asmbPLS.cv <- function(X.matrix, Y.matrix, PLS.comp, X.dim, quantile.comb.table, Y.indicator, k = 5, only.observe = TRUE, seed = 1) {
  ## error check
  stopifnot(!missing(X.matrix),
            !missing(Y.matrix),
            !missing(PLS.comp),
            !missing(X.dim),
            !missing(quantile.comb.table),
            !missing(Y.indicator),
            is.matrix(X.matrix),
            is.matrix(Y.matrix),
            is.matrix(quantile.comb.table),
            is.numeric(PLS.comp))
  return(asmbPLS_CV(X.matrix, Y.matrix, PLS.comp, X.dim, quantile.comb.table, Y.indicator, k, only.observe, seed))
}
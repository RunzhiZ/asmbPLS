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
#' @param ncv The number of repetitions of CV. The default is 5.
#' @param only.observe Whether only observed samples in the validation set 
#' should be used for calculating the MSE for CV. The default is TRUE.
#' @param expected.measure.decrease The measure you expect to decrease by percent 
#' after including one more PLS component, which will affect the selection of optimal 
#' PLS components. The default is 0.05 (5\%).
#' @param center A logical value indicating whether mean center should be 
#' implemented for X.matrix and Y.matrix. The default is TRUE.
#' @param scale  A logical value indicating whether scale should be 
#' implemented for X.matrix and Y.matrix. The default is TRUE.
#' @param maxiter A integer indicating the maximum number of iteration. The
#' default number is 100.
#' 
#' @return 
#' \code{asmbPLS.cv} returns a list containing the following components:
#' \item{quantile_table_CV}{A matrix containing the selected quantile 
#' combination and the corresponding measures of CV for each PLS component.}
#' \item{optimal_nPLS}{Optimal number of PLS components.}.
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
#'                          ncv = 3)
#' quantile.comb <- cv.results$quantile_table_CV[,1:length(X.dim)]
#' n.PLS <- cv.results$optimal_nPLS
#'  
#' ## asmbPLS fit
#' asmbPLS.results <- asmbPLS.fit(X.matrix = X.matrix, 
#'                                Y.matrix = Y.matrix, 
#'                                PLS.comp = n.PLS, 
#'                                X.dim = X.dim, 
#'                                quantile.comb = quantile.comb)
#' 
#' @export
#' @useDynLib asmbPLS, .registration=TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom stats quantile

asmbPLS.cv <- function(X.matrix, Y.matrix, PLS.comp, X.dim, quantile.comb.table, 
                       Y.indicator, k = 5, ncv = 5, only.observe = TRUE, 
                       expected.measure.decrease = 0.05, 
                       center = TRUE, scale = TRUE, maxiter = 100) {
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
  cv_results <- asmbPLS_CV(X.matrix, Y.matrix, PLS.comp, X.dim, quantile.comb.table, Y.indicator, k, ncv, only.observe, expected.measure.decrease, center, scale, maxiter)
  colnames(cv_results$quantile_table_CV) <- c(paste0("X.", 1:length(X.dim)), "MSE")
  return(cv_results)
}
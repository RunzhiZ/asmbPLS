#' Using an asmbPLS model for prediction of new samples
#'
#' Derives predictions for new samples from a model fitted by the function
#' \code{\link[asmbPLS]{asmbPLS.fit}} or \code{\link[asmbPLS]{mbPLS.fit}}.
#' 
#' @param fit.results The output of either \code{\link[asmbPLS]{asmbPLS.fit}} or 
#' \code{\link[asmbPLS]{mbPLS.fit}}.
#' @param X.matrix.new A predictors matrix, whose predictors are the same as 
#' the predictors in model fitting.
#' @param PLS.comp Number of PLS components used for prediction.
#' 
#' @return 
#' \code{asmbPLSDA.predict} returns a list containing the following components:
#' \item{Y_pred}{Predicted value for the new sampels.}
#' \item{NewX_super_score}{Predicted super score for new samples, which can be
#' used as predictors for other regression models.}
#' 
#' @examples
#' ## Use the example dataset
#' data(asmbPLS.example)
#' X.matrix = asmbPLS.example$X.matrix
#' X.matrix.new = asmbPLS.example$X.matrix.new
#' Y.matrix = asmbPLS.example$Y.matrix
#' PLS.comp = asmbPLS.example$PLS.comp
#' X.dim = asmbPLS.example$X.dim
#' quantile.comb = asmbPLS.example$quantile.comb
#'  
#' ## asmbPLS fit
#' asmbPLS.results <- asmbPLS.fit(X.matrix = X.matrix, 
#'                                Y.matrix = Y.matrix, 
#'                                PLS.comp = PLS.comp, 
#'                                X.dim = X.dim, 
#'                                quantile.comb = quantile.comb)
#' 
#' ## asmbPLS prediction for the new data, you could use different numbers of 
#' ## PLS components for prediction
#' ## Use only the first PLS component 
#' Y.pred.1 <- asmbPLS.predict(asmbPLS.results, X.matrix.new, 1)
#' ## Use the first two PLS components
#' Y.pred.2 <- asmbPLS.predict(asmbPLS.results, X.matrix.new, 2)
#' 
#' @export
#' @useDynLib asmbPLS, .registration=TRUE
#' @importFrom Rcpp sourceCpp

asmbPLS.predict <- function(fit.results, X.matrix.new, PLS.comp){
  stopifnot(!missing(X.matrix.new),
            !missing(PLS.comp),
            !missing(fit.results),
            is.matrix(X.matrix.new))
  return(asmbPLS_predict(fit.results, X.matrix.new, PLS.comp))
}
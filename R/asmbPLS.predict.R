#' Using an asmbPLS model for prediction of new observations
#'
#' Derives predictions for new observations from a model fitted by the function
#' \code{\link[asmbPLS]{asmbPLS.fit}} or \code{\link[asmbPLS]{mbPLS.fit}}.
#' 
#' @param X.matrix.new A predictors matrix, whose predictors are the same as 
#' the predictors in model fitting.
#' @param PLS.comp Number of PLS components used for prediction.
#' @param fit.results The output of either \code{\link[asmbPLS]{asmbPLS.fit}} or 
#' \code{\link[asmbPLS]{mbPLS.fit}}.
#' 
#' @return 
#' \code{asmbPLS.predict} returns a matrix containing the prediction for the 
#' new data.
#' 
#' @examples
#' ## Use the example dataset
#' data(asmbPLS.predict.example)
#'  
#' ## asmbPLS fit
#' asmbPLS.results <- asmbPLS.fit(X.matrix = X.matrix, Y.matrix = Y.matrix, 
#' PLS.comp = PLS.comp, X.dim = X.dim, quantile.comb = quantile.comb)
#' 
#' ## asmbPLS prediction for the new data, you could use different numbers of 
#' ## PLS components for prediction
#' Y.pred.1 <- asmbPLS.predict(X.matrix.new, 1, asmbPLS.results)
#' Y.pred.2 <- asmbPLS.predict(X.matrix.new, 2, asmbPLS.results)
#' 
#' @export
#' @useDynLib asmbPLS, .registration=TRUE
#' @importFrom Rcpp sourceCpp

asmbPLS.predict <- function(X.matrix.new, PLS.comp, fit.results){
  stopifnot(!missing(X.matrix.new),
            !missing(PLS.comp),
            !missing(fit.results),
            is.matrix(X.matrix.new))
  return(asmbPLS_predict_rcpp(X.matrix.new, PLS.comp, fit.results))
}
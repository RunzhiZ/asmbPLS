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
#' \code{asmbPLS.predict} returns a matrix containing the prediction for the 
#' new data.
#' 
#' @examples
#' ## Use the example dataset
#' data(asmbPLS.predict.example)
#'  
#' ## asmbPLS fit
#' asmbPLS.results <- asmbPLS.fit(
#' X.matrix = asmbPLS.predict.example$X.matrix, 
#' Y.matrix = asmbPLS.predict.example$Y.matrix, 
#' PLS.comp = asmbPLS.predict.example$PLS.comp, 
#' X.dim = asmbPLS.predict.example$X.dim, 
#' quantile.comb = asmbPLS.predict.example$quantile.comb)
#' 
#' ## asmbPLS prediction for the new data, you could use different numbers of 
#' ## PLS components for prediction
#' Y.pred.1 <- asmbPLS.predict(asmbPLS.results, asmbPLS.predict.example$X.matrix.new, 1)
#' Y.pred.2 <- asmbPLS.predict(asmbPLS.results, asmbPLS.predict.example$X.matrix.new, 2)
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
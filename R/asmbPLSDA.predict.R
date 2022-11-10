#' Using an asmbPLS-DA model for prediction of new samples
#'
#' Derives predictions for new samples from a model fitted by the function
#' \code{\link[asmbPLS]{asmbPLSDA.fit}}.
#' 
#' @param X.matrix.new A predictors matrix, whose predictors are the same as 
#' the predictors in model fitting.
#' @param PLS.comp Number of PLS components used for prediction.
#' @param fit.results The output of either \code{\link[asmbPLS]{asmbPLSDA.fit}}
#' @param Method Decision rule used for prediction. For binary outcome, the 
#' methods include "\code{fixed_cutoff}" (default), "\code{Euclidean_distance_X}" 
#' and "\code{Mahalanobis_distance_X}". For categorical outcome with more than 2 
#' levels, the methods include "\code{Max_Y}" (default), 
#' "\code{Euclidean_distance_X}", "\code{Mahalanobis_distance_X}", 
#' "\code{Euclidean_distance_Y}", and "\code{PCA_Mahalanobis_distance_Y}".
#' 
#' @return 
#' \code{asmbPLSDA.predict} returns a list containing the following components:
#' \item{Y_pred}{Predicted class for the new sampels.}
#' \item{Y_pred_numeric}{Predicted Y values for the new samples, different 
#' decision rules can be used to obtain different Y_pred.}
#' \item{NewX_super_score}{Predicted super score for new samples, which can be
#' used as predictors for other classification algorithms.}
#' \item{Method}{Decision rule used for preidction.}
#' 
#' @examples
#' ## Use the example dataset
#' data(asmbPLSDA.predict.example)
#'  
#' ## asmbPLSDA fit for binary outcome
#' asmbPLSDA.binary.results <- asmbPLSDA.fit(
#' X.matrix = asmbPLSDA.predict.example$X.matrix, 
#' Y.matrix = asmbPLSDA.predict.example$Y.matrix.binary, 
#' PLS.comp = asmbPLSDA.predict.example$PLS.comp, 
#' X.dim = asmbPLSDA.predict.example$X.dim, 
#' quantile.comb = asmbPLSDA.predict.example$quantile.comb,
#' outcome.type = "binary")
#' 
#' ## asmbPLSDA fit for categorical outcome with more than 2 levels
#' asmbPLSDA.morethan2levels.results <- asmbPLSDA.fit(
#' X.matrix = asmbPLSDA.predict.example$X.matrix, 
#' Y.matrix = asmbPLSDA.predict.example$Y.matrix.morethan2levels, 
#' PLS.comp = asmbPLSDA.predict.example$PLS.comp, 
#' X.dim = asmbPLSDA.predict.example$X.dim, 
#' quantile.comb = asmbPLSDA.predict.example$quantile.comb,
#' outcome.type = "morethan2levels")
#' 
#' ## asmbPLSDA prediction for the new data, you could use different numbers of 
#' ## PLS components for prediction
#' Y.pred.binary.1 <- asmbPLSDA.predict(
#' asmbPLSDA.predict.example$X.matrix.new, 1, 
#' asmbPLSDA.binary.results)
#' Y.pred.binary.2 <- asmbPLSDA.predict(
#' asmbPLSDA.predict.example$X.matrix.new, 2, 
#' asmbPLSDA.binary.results)
#' Y.pred.morethan2levels.1 <- asmbPLSDA.predict(
#' asmbPLSDA.predict.example$X.matrix.new, 1, 
#' asmbPLSDA.morethan2levels.results)
#' Y.pred.morethan2levels.2 <- asmbPLSDA.predict(
#' asmbPLSDA.predict.example$X.matrix.new, 2, 
#' asmbPLSDA.morethan2levels.results)
#' 
#' @export
#' @useDynLib asmbPLS, .registration=TRUE
#' @importFrom Rcpp sourceCpp

asmbPLSDA.predict <- function(X.matrix.new, PLS.comp, fit.results, Method = NULL){
  stopifnot(!missing(X.matrix.new),
            !missing(PLS.comp),
            !missing(fit.results),
            is.matrix(X.matrix.new))
  return(asmbPLSDA_predict(X.matrix.new, PLS.comp, fit.results, Method))
}
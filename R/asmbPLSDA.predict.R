#' Using an asmbPLS-DA model for classification of new samples
#'
#' Derives classification for new samples from a model fitted by the function
#' \code{\link[asmbPLS]{asmbPLSDA.fit}}.
#' 
#' @param fit.results The output of \code{\link[asmbPLS]{asmbPLSDA.fit}}
#' @param X.matrix.new A predictors matrix, whose predictors are the same as 
#' the predictors in model fitting.
#' @param PLS.comp Number of PLS components used for prediction.
#' @param method Decision rule used for prediction. For binary outcome, the 
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
#' \item{method}{Decision rule used for preidction.}
#' 
#' @examples
#' ## Use the example dataset
#' data(asmbPLSDA.example)
#' X.matrix = asmbPLSDA.example$X.matrix
#' X.matrix.new = asmbPLSDA.example$X.matrix.new
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
#' ## asmbPLSDA prediction for the new data, you could use different numbers of 
#' ## PLS components for prediction
#' ## Use only the first PLS component 
#' Y.pred.binary.1 <- asmbPLSDA.predict(asmbPLSDA.fit.binary, 
#'                                      X.matrix.new, 
#'                                      PLS.comp = 1)
#' ## Use the first two PLS components                                      
#' Y.pred.binary.2 <- asmbPLSDA.predict(asmbPLSDA.fit.binary,
#'                                      X.matrix.new, 
#'                                      PLS.comp = 2)
#' 
#' ## PLS components for prediction
#' Y.pred.multiclass.1 <- asmbPLSDA.predict(asmbPLSDA.fit.multiclass,
#'                                          X.matrix.new, 
#'                                          PLS.comp = 1)
#' ## Use the first two PLS components     
#' Y.pred.multiclass.2 <- asmbPLSDA.predict(asmbPLSDA.fit.multiclass,
#'                                          X.matrix.new, 
#'                                          PLS.comp = 2)
#' 
#' @export
#' @useDynLib asmbPLS, .registration=TRUE
#' @importFrom Rcpp sourceCpp

asmbPLSDA.predict <- function(fit.results, X.matrix.new, PLS.comp, method = NULL){
  outcome.type <- fit.results$Outcome_type
  if(outcome.type == "binary" & is.null(method)) {method <- "fixed_cutoff"}
  if(outcome.type == "multiclass" & is.null(method)) {method <- "Max_Y"}
  stopifnot(!missing(X.matrix.new),
            !missing(PLS.comp),
            !missing(fit.results),
            is.matrix(X.matrix.new))
  return(asmbPLSDA_predict(fit.results, X.matrix.new, PLS.comp, method))
}
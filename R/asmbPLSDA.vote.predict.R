#' Using an asmbPLS-DA vote model for classification of new samples
#'
#' Function to make the classification using the weights and fitted model 
#' obtained from \code{\link{asmbPLSDA.vote.fit}}. The final classification 
#' results are the weighted classification using the decision rules included.
#' 
#' @param fit.results The output of \code{\link{asmbPLSDA.vote.fit}}.
#' @param X.matrix.new A predictors matrix, whose predictors are the same as 
#' the predictors in model fitting.
#' 
#' @return 
#' \item{Y_pred}{Predicted class for the new sampels.}

#' 
#' @examples
#' ## Use the example dataset
#' data(asmbPLSDA.example)
#' X.matrix = asmbPLSDA.example$X.matrix
#' X.matrix.new = asmbPLSDA.example$X.matrix.new
#' Y.matrix.binary = asmbPLSDA.example$Y.matrix.binary
#' X.dim = asmbPLSDA.example$X.dim
#' PLS.comp = asmbPLSDA.example$PLS.comp
#' quantile.comb.table.cv = asmbPLSDA.example$quantile.comb.table.cv
#' 
#' ## Cross validaiton based on fixed cutoff
#' cv.results.cutoff <- asmbPLSDA.cv(X.matrix = X.matrix, 
#'                                   Y.matrix = Y.matrix.binary,
#'                                   PLS.comp = PLS.comp, 
#'                                   X.dim = X.dim, 
#'                                   quantile.comb.table = quantile.comb.table.cv, 
#'                                   outcome.type = "binary", 
#'                                   method = "fixed_cutoff",
#'                                   k = 3,
#'                                   ncv = 1)
#' quantile.comb.cutoff <- cv.results.cutoff$quantile_table_CV
#' 
#' ## Cross validation using Euclidean distance of X super score
#' cv.results.EDX <- asmbPLSDA.cv(X.matrix = X.matrix, 
#'                                Y.matrix = Y.matrix.binary,
#'                                PLS.comp = PLS.comp, 
#'                                X.dim = X.dim, 
#'                                quantile.comb.table = quantile.comb.table.cv, 
#'                                outcome.type = "binary", 
#'                                method = "Euclidean_distance_X",
#'                                k = 3,
#'                                ncv = 1)
#' quantile.comb.EDX <- cv.results.EDX$quantile_table_CV
#' 
#' ## Cross validation using Mahalanobis distance of X super score
#' cv.results.MDX <- asmbPLSDA.cv(X.matrix = X.matrix, 
#'                                   Y.matrix = Y.matrix.binary,
#'                                   PLS.comp = PLS.comp, 
#'                                   X.dim = X.dim, 
#'                                   quantile.comb.table = quantile.comb.table.cv, 
#'                                   outcome.type = "binary", 
#'                                   method = "Mahalanobis_distance_X",
#'                                   k = 3,
#'                                   ncv = 1)
#' quantile.comb.MDX <- cv.results.MDX$quantile_table_CV
#' 
#' #### vote list ####
#' cv.results.list = list(fixed_cutoff = quantile.comb.cutoff,
#'                        Euclidean_distance_X = quantile.comb.EDX,
#'                        Mahalanobis_distance_X = quantile.comb.MDX)
#' 
#' ## vote models fit
#' vote.fit <- asmbPLSDA.vote.fit(X.matrix = X.matrix, 
#'                                Y.matrix = Y.matrix.binary, 
#'                                X.dim = X.dim, 
#'                                nPLS = c(cv.results.cutoff$optimal_nPLS, 
#'                                cv.results.EDX$optimal_nPLS, 
#'                                cv.results.MDX$optimal_nPLS),
#'                                cv.results.list = cv.results.list, 
#'                                outcome.type = "binary",
#'                                method = "weighted")
#' 
#' ## classification
#' vote.predict <- asmbPLSDA.vote.predict(vote.fit, X.matrix.new)
#' 
#' @export
#' @useDynLib asmbPLS, .registration=TRUE
#' @importFrom Rcpp sourceCpp

asmbPLSDA.vote.predict <- function(fit.results, X.matrix.new) {
  G <- ncol(fit.results[[1]]$fit.model$Y_group)
  outcome.type <- fit.results[[1]]$outcome.type
  predict.results <- matrix(0, nrow = nrow(X.matrix.new), ncol = G)
  predict.results.output <- matrix(0, nrow = nrow(X.matrix.new), ncol = G)
  for(i in 1:length(fit.results)) {
    fit.single <- fit.results[[i]]
    method_name = names(fit.results)[i]
    predict.single <- asmbPLSDA.predict(fit.single$fit.model, X.matrix.new, fit.single$nPLS, method = method_name)
    predict.results = predict.results + fit.single$weight * predict.single$Y_pred
  }
  
  if(outcome.type == "binary"){
    predict.results.output <- matrix(as.numeric(predict.results > 0.5))
  }
  if(outcome.type == "multiclass"){
    for(i in 1:nrow(predict.results)) {
      max_index <- which.max(predict.results[i,])
      predict.results.output[i, max_index] <- 1
    }
  }
  return(predict.results.output)
}
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
#' Y.matrix.morethan2levels = asmbPLSDA.example$Y.matrix.morethan2levels
#' X.dim = asmbPLSDA.example$X.dim
#' PLS.comp = asmbPLSDA.example$PLS.comp
#' quantile.comb.table.cv = asmbPLSDA.example$quantile.comb.table.cv
#' 
#' ## Cross validaiton based on max Y
#' cv.results.max <- asmbPLSDA.cv(X.matrix = X.matrix, 
#'                                Y.matrix = Y.matrix.morethan2levels,
#'                                PLS.comp = 5, 
#'                                X.dim = X.dim, 
#'                                quantile.comb.table = quantile.comb.table.cv, 
#'                                outcome.type = "morethan2levels", 
#'                                Method = "Max_Y")
#' quantile.comb.max <- cv.results.max$quantile_table_CV
#' 
#' ## Cross validation using Euclidean distance of X super score
#' cv.results.EDX <- asmbPLSDA.cv(X.matrix = X.matrix, 
#'                                Y.matrix = Y.matrix.morethan2levels,
#'                                PLS.comp = 5, 
#'                                X.dim = X.dim, 
#'                                quantile.comb.table = quantile.comb.table.cv, 
#'                                outcome.type = "morethan2levels", 
#'                                Method = "Euclidean_distance_X")
#' quantile.comb.EDX <- cv.results.EDX$quantile_table_CV
#' 
#' ## Cross validation using PCA + Mahalanobis distance of Y
#' cv.results.PCAMDY <- asmbPLSDA.cv(X.matrix = X.matrix, 
#'                                   Y.matrix = Y.matrix.morethan2levels,
#'                                   PLS.comp = 5, 
#'                                   X.dim = X.dim, 
#'                                   quantile.comb.table = quantile.comb.table.cv, 
#'                                   outcome.type = "morethan2levels", 
#'                                   Method = "Euclidean_distance_X")
#' quantile.comb.PCAMDY <- cv.results.PCAMDY$quantile_table_CV
#' 
#' #### vote list ####
#' cv.results.list = list(Max_Y = quantile.comb.max,
#'                        Euclidean_distance_X = quantile.comb.EDX,
#'                        PCA_Mahalanobis_distance_Y = quantile.comb.PCAMDY)
#' 
#' ## vote models fit
#' vote.fit <- asmbPLSDA.vote.fit(X.matrix = X.matrix, 
#'                                Y.matrix = Y.matrix.morethan2levels, 
#'                                X.dim = X.dim, 
#'                                nPLS = c(2, 2, 2),
#'                                cv.results.list = cv.results.list, 
#'                                outcome.type = "morethan2levels",
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
    Method_name = names(fit.results)[i]
    predict.single <- asmbPLSDA.predict(fit.single$fit.model, X.matrix.new, fit.single$nPLS, Method = Method_name)
    predict.results = predict.results + fit.single$weight * predict.single$Y_pred
  }
  
  if(outcome.type == "binary"){
    predict.results.output <- matrix(as.numeric(predict.results > 0.5))
  }
  if(outcome.type == "morethan2levels"){
    for(i in 1:nrow(predict.results)) {
      max_index <- which.max(predict.results[i,])
      predict.results.output[i, max_index] <- 1
    }
  }
  return(predict.results.output)
}
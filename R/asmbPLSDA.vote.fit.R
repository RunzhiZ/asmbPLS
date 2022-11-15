#' asmbPLS-DA vote model fit 
#'
#' Function to fit multiple asmbPLS-DA models using cross validation results with 
#' different decision rules obtained from \code{\link[asmbPLS]{asmbPLSDA.cv}}, 
#' the weight for each model are calculated based on cross-validation accuracy, 
#' which can be used for \code{\link[asmbPLS]{asmbPLSDA.vote.predict}} to obtain
#' the final classification.
#' 
#' @param X.matrix Predictors matrix. Samples in rows, variables in columns.
#' @param Y.matrix Outcome matrix. Samples in rows, this is a matrix with one 
#' column (binary) or multiple columns (more than 2 levels, dummy variables).
#' @param X.dim A vector containing the number of predictors in each block 
#' (ordered).
#' @param cv.results.list A list containing \code{quantile_table_CV} from 
#' \code{\link[asmbPLS]{asmbPLSDA.cv}} using different decision rules, the name
#' of each element in the list should be the corresponding name of decision 
#' rule. 
#' @param outcome.type The type of the outcome Y. "\code{binary}" for binary 
#' outcome, and "\code{morethan2levels}" for categorical outcome with more than 2 
#' levels.
#' @param expected.accuracy.increase The accuracy you expect to increase after 
#' including one more PLS component, which will affect the selection of optimal 
#' PLS components. The default is 0.005.
#' @param center A logical value indicating whether weighted mean center should 
#' be implemented for \code{X.matrix} and \code{Y.matrix}. The default is TRUE.
#' @param scale  A logical value indicating whether scale should be 
#' implemented for \code{X.matrix}. The default is TRUE.
#' 
#' @return 
#' \code{asmbPLSDA.vote.fit} returns a list of lists, which can be used as the 
#' inputs for \code{\link{asmbPLSDA.vote.predict}}. Each list contains the 
#' fit information for model with specific decision rule: 
#' \item{fit.model}{A list containing model fit information.}
#' \item{nPLS.optimal}{The optimal number of PLS components.}
#' \item{weight}{The weight for this model.}
#' 
#' @examples
#' ## Use the example dataset
#' data(asmbPLSDA.example)
#' X.matrix = asmbPLSDA.example$X.matrix
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
#'                                cv.results.list = cv.results.list, 
#'                                outcome.type = "morethan2levels")
#' 
#' @export
#' @useDynLib asmbPLS, .registration=TRUE
#' @importFrom Rcpp sourceCpp

asmbPLSDA.vote.fit <- function(X.matrix, 
                               Y.matrix, 
                               X.dim, 
                               cv.results.list, 
                               outcome.type,
                               expected.accuracy.increase = 0.005,
                               center = TRUE, 
                               scale = TRUE) {
  
  n_dim = length(X.dim)
  fit.list = cv.results.list
  accuracy.weight <- NULL
  nPLS_optimal <- NULL
  
  for(i in 1:length(cv.results.list)) {
    Method_used = names(cv.results.list)[i]
    cv.results <- cv.results.list[[i]]
    accuracy = cv.results[, n_dim + 1]
    ## obtain optimal nPLS
    nPLS_optimal_temp = 1
    if(length(accuracy) > 1) {
      for(j in 1:(length(accuracy) - 1)) {
        if(accuracy[j + 1] > (accuracy[j] + expected.accuracy.increase)) {
          nPLS_optimal_temp = nPLS_optimal_temp + 1
        }
      }
    }
    nPLS_optimal[i] <- nPLS_optimal_temp
    accuracy.weight[i] <- accuracy[nPLS_optimal_temp]
  }
  ## obtain weight for each method
  accuracy.weight = log(accuracy.weight/(1-accuracy.weight))
  accuracy.weight = accuracy.weight/sum(accuracy.weight)
  
  for(i in 1:length(cv.results.list)) {
    cv.results <- cv.results.list[[i]]
    quantile.comb = matrix(cv.results[1:nPLS_optimal[i], 1:n_dim], nrow = nPLS_optimal[i])
    fit.model <- asmbPLSDA.fit(X.matrix,
                               Y.matrix,
                               nPLS_optimal[i],
                               X.dim,
                               quantile.comb,
                               outcome.type,
                               center,
                               scale)
    fit.list[[i]] <- list(fit.model = fit.model, 
                          nPLS.optimal = nPLS_optimal[i], 
                          weight = accuracy.weight[i])
  }
  return(fit.list)
}
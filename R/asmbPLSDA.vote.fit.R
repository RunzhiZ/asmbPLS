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
#' @param nPLS A vector containing the number of PLS components used for 
#' different decision rules.
#' @param outcome.type The type of the outcome Y. "\code{binary}" for binary 
#' outcome, and "\code{multiclass}" for categorical outcome with more than 2 
#' levels.
#' @param method Vote options. "\code{unweighted}" gives each decision rule the 
#' same weight; "\code{weighted}" assigns higher weight to method with higher 
#' measure, i.e. weight = log(measure/(1-measure)); "\code{ranked}" ranks
#' the given methods based on the average rank of methods using accuracy, 
#' balanced accuracy, precision, recall and F1 score.
#' @param measure Measure to be selected when \code{method} is \code{weighted}. 
#' The default is \code{B_accuracy}.
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
#' \item{nPLS}{The number of PLS components used.}
#' \item{weight}{The weight for this model.}
#' \item{outcome.type}{The type of the outcome Y.}
#' 
#' @examples
#' ## Use the example dataset
#' data(asmbPLSDA.example)
#' X.matrix = asmbPLSDA.example$X.matrix
#' X.matrix.new = asmbPLSDA.example$X.matrix.new
#' Y.matrix.multiclass = asmbPLSDA.example$Y.matrix.morethan2levels
#' X.dim = asmbPLSDA.example$X.dim
#' PLS.comp = asmbPLSDA.example$PLS.comp
#' quantile.comb.table.cv = asmbPLSDA.example$quantile.comb.table.cv
#' 
#' ## Cross validaiton based on max Y
#' cv.results.max <- asmbPLSDA.cv(X.matrix = X.matrix, 
#'                                Y.matrix = Y.matrix.multiclass,
#'                                PLS.comp = 3, 
#'                                X.dim = X.dim, 
#'                                quantile.comb.table = quantile.comb.table.cv, 
#'                                outcome.type = "multiclass", 
#'                                method = "Max_Y")
#' quantile.comb.max <- cv.results.max$quantile_table_CV
#' 
#' ## Cross validation using Euclidean distance of X super score
#' cv.results.EDX <- asmbPLSDA.cv(X.matrix = X.matrix, 
#'                                Y.matrix = Y.matrix.multiclass,
#'                                PLS.comp = 3, 
#'                                X.dim = X.dim, 
#'                                quantile.comb.table = quantile.comb.table.cv, 
#'                                outcome.type = "multiclass", 
#'                                method = "Euclidean_distance_X")
#' quantile.comb.EDX <- cv.results.EDX$quantile_table_CV
#' 
#' ## Cross validation using PCA + Mahalanobis distance of Y
#' cv.results.PCAMDY <- asmbPLSDA.cv(X.matrix = X.matrix, 
#'                                   Y.matrix = Y.matrix.multiclass,
#'                                   PLS.comp = 3, 
#'                                   X.dim = X.dim, 
#'                                   quantile.comb.table = quantile.comb.table.cv, 
#'                                   outcome.type = "multiclass", 
#'                                   method = "Euclidean_distance_X")
#' quantile.comb.PCAMDY <- cv.results.PCAMDY$quantile_table_CV
#' 
#' #### vote list ####
#' cv.results.list = list(Max_Y = quantile.comb.max,
#'                        Euclidean_distance_X = quantile.comb.EDX,
#'                        PCA_Mahalanobis_distance_Y = quantile.comb.PCAMDY)
#' 
#' ## vote models fit
#' vote.fit <- asmbPLSDA.vote.fit(X.matrix = X.matrix, 
#'                                Y.matrix = Y.matrix.multiclass, 
#'                                X.dim = X.dim, 
#'                                nPLS = c(cv.results.max$optimal_nPLS, 
#'                                cv.results.EDX$optimal_nPLS, 
#'                                cv.results.PCAMDY$optimal_nPLS),
#'                                cv.results.list = cv.results.list, 
#'                                outcome.type = "multiclass",
#'                                method = "weighted")
#' 
#' @export
#' @useDynLib asmbPLS, .registration=TRUE
#' @importFrom Rcpp sourceCpp

asmbPLSDA.vote.fit <- function(X.matrix, 
                               Y.matrix, 
                               X.dim, 
                               cv.results.list,
                               nPLS,
                               outcome.type,
                               method = "weighted",
                               measure = NULL,
                               center = TRUE, 
                               scale = TRUE) {
  
  n_dim = length(X.dim)
  fit.list = cv.results.list
  weight <- rep(0, length(cv.results.list))
  measure <- NULL
  
  
  if (method == "weighted") {
    if (is.null(measure)) {measure <- "B_accuracy"}
    measure_list <- c("accuracy", "B_accuracy", "precision", "recall", "F1")
    measure_indices <- 1:5
    measure_index <- measure_indices[which(measure_list == measure)]
    for(i in 1:length(cv.results.list)) {
      cv.results <- cv.results.list[[i]]
      measure_value = cv.results[, n_dim + 1]
      weight[i] <- measure_value[nPLS[i]]
    }
    ## obtain weight for each method
    weight = log(weight/(1-weight))
    weight = weight/sum(weight)
  }
  if (method == "unweighted") {
    for(i in 1:length(cv.results.list)) {
      weight[i] <- 1/length(cv.results.list)
    }
  }
  if (method == "ranked") {
    for(i in 1:length(cv.results.list)) {
      cv.results <- cv.results.list[[i]]
      measure <- rbind(measure, cv.results[nPLS[i], (n_dim + 1):(n_dim + 4)])
    }
    rank_average <- apply(apply(measure, 2, rank), 1, mean)
    max_index <- which.max(rank_average)
    weight[max_index] <- 1
  }
  
  for(i in 1:length(cv.results.list)) {
    cv.results <- cv.results.list[[i]]
    quantile.comb = matrix(cv.results[1:nPLS[i], 1:n_dim], nrow = nPLS[i])
    fit.model <- asmbPLSDA.fit(X.matrix,
                               Y.matrix,
                               nPLS[i],
                               X.dim,
                               quantile.comb,
                               outcome.type,
                               center,
                               scale)
    fit.list[[i]] <- list(fit.model = fit.model, 
                          nPLS = nPLS[i], 
                          weight = weight[i],
                          outcome.type = outcome.type)
  }
  
  return(fit.list)
}
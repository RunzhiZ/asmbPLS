#' Cross-validation for asmbPLS-DA to find the best combinations of quantiles for classification
#'
#' Function to find the best combinations of quantiles used for classification via
#' cross-validation. Usually should be conducted before 
#' \code{\link[asmbPLS]{asmbPLSDA.fit}} to obtain the quantile combinations.
#' 
#' @param X.matrix Predictors matrix. Samples in rows, variables in columns.
#' @param Y.matrix Outcome matrix. Samples in rows, this is a matrix with one 
#' column (binary) or multiple columns (more than 2 levels, dummy variables).
#' @param PLS.comp Number of PLS components in asmbPLS-DA.
#' @param X.dim A vector containing the number of predictors in each block 
#' (ordered).
#' @param quantile.comb.table A matrix containing user-defined quantile 
#' combinations used for CV, whose column number equals to the 
#' number of blocks.
#' @param outcome.type The type of the outcome Y. "\code{binary}" for binary 
#' outcome, and "\code{multiclass}" for categorical outcome with more than 2 
#' levels.
#' @param method Decision rule used for CV. For binary outcome, the 
#' methods include "\code{fixed_cutoff}", "\code{Euclidean_distance_X}" and
#' "\code{Mahalanobis_distance_X}". For categorical outcome with more than 2 
#' levels, the methods include "\code{Max_Y}", "\code{Euclidean_distance_X}",
#' "\code{Mahalanobis_distance_X}", "\code{Euclidean_distance_Y}", and 
#' "\code{PCA_Mahalanobis_distance_Y}".
#' @param measure Five measures are available: overall accuracy \code{accuracy},
#' balanced accuracy \code{B_accuracy}, precision \code{precision}, recall
#' \code{recall}, F1 score \code{F1}.
#' @param k The number of folds of CV procedure. The default is 5.
#' @param ncv The number of repetitions of CV. The default is 5.
#' @param expected.measure.increase The measure you expect to increase after 
#' including one more PLS component, which will affect the selection of optimal 
#' PLS components. The default is 0.005.
#' @param center A logical value indicating whether weighted mean center should 
#' be implemented for \code{X.matrix} and \code{Y.matrix}. The default is TRUE.
#' @param scale  A logical value indicating whether scale should be 
#' implemented for \code{X.matrix}. The default is TRUE.
#' @param maxiter A integer indicating the maximum number of iteration. The
#' default number is 100.
#' 
#' @return 
#' \code{asmbPLSDA.cv} returns a list containing the following components:
#' \item{quantile_table_CV}{A matrix containing the selected quantile 
#' combination and the corresponding measures of CV for each PLS component.}
#' \item{optimal_nPLS}{Optimal number of PLS components.}
#' 
#' @examples
#' ## Use the example dataset
#' data(asmbPLSDA.example)
#' X.matrix = asmbPLSDA.example$X.matrix
#' Y.matrix.binary = asmbPLSDA.example$Y.matrix.binary
#' Y.matrix.multiclass = asmbPLSDA.example$Y.matrix.morethan2levels
#' X.dim = asmbPLSDA.example$X.dim
#' PLS.comp = asmbPLSDA.example$PLS.comp
#' quantile.comb.table.cv = asmbPLSDA.example$quantile.comb.table.cv
#' 
#' ## cv to find the best quantile combinations for model fitting (binary outcome)
#' cv.results.binary <- asmbPLSDA.cv(X.matrix = X.matrix, 
#'                                   Y.matrix = Y.matrix.binary, 
#'                                   PLS.comp = PLS.comp, 
#'                                   X.dim = X.dim, 
#'                                   quantile.comb.table = quantile.comb.table.cv, 
#'                                   outcome.type = "binary")
#' quantile.comb.binary <- cv.results.binary$quantile_table_CV[,1:length(X.dim)]
#' 
#' ## asmbPLSDA fit using the selected quantile combination (binary outcome)
#' asmbPLSDA.fit.binary <- asmbPLSDA.fit(X.matrix = X.matrix, 
#'                                       Y.matrix = Y.matrix.binary, 
#'                                       PLS.comp = PLS.comp, 
#'                                       X.dim = X.dim, 
#'                                       quantile.comb = quantile.comb.binary,
#'                                       outcome.type = "binary")
#' 
#' 
#' ## cv to find the best quantile combinations for model fitting 
#' ## (categorical outcome with more than 2 levels)
#' cv.results.multiclass <- asmbPLSDA.cv(X.matrix = X.matrix, 
#'                                       Y.matrix = Y.matrix.multiclass, 
#'                                       PLS.comp = PLS.comp, 
#'                                       X.dim = X.dim, 
#'                                       quantile.comb.table = quantile.comb.table.cv, 
#'                                       outcome.type = "multiclass")
#' quantile.comb.multiclass <- cv.results.multiclass$quantile_table_CV[,1:length(X.dim)]
#' 
#' ## asmbPLSDA fit (categorical outcome with more than 2 levels)
#' asmbPLSDA.fit.multiclass <- asmbPLSDA.fit(X.matrix = X.matrix, 
#'                                           Y.matrix = Y.matrix.multiclass, 
#'                                           PLS.comp = PLS.comp, 
#'                                           X.dim = X.dim, 
#'                                           quantile.comb = quantile.comb.multiclass,
#'                                           outcome.type = "multiclass")
#' 
#' @export
#' @useDynLib asmbPLS, .registration=TRUE
#' @importFrom Rcpp sourceCpp

asmbPLSDA.cv <- function(X.matrix, Y.matrix, PLS.comp, X.dim, 
                         quantile.comb.table, outcome.type, method = NULL, 
                         measure = "B_accuracy", k = 5, ncv = 5, 
                         expected.measure.increase = 0.005, 
                         center = TRUE, scale = TRUE, maxiter = 100) {
  if(outcome.type == "binary" & is.null(method)) {method <- "fixed_cutoff"}
  if(outcome.type == "multiclass" & is.null(method)) {method <- "Max_Y"}
  ## error check
  stopifnot(!missing(X.matrix),
            !missing(Y.matrix),
            !missing(PLS.comp),
            !missing(X.dim),
            !missing(quantile.comb.table),
            !missing(outcome.type),
            is.matrix(X.matrix),
            is.matrix(Y.matrix),
            is.matrix(quantile.comb.table),
            is.numeric(PLS.comp))
  cv_results <- asmbPLSDA_CV(X.matrix, Y.matrix, PLS.comp, X.dim, quantile.comb.table, outcome.type, method, measure, k, ncv, expected.measure.increase, center, scale, maxiter)
  colnames(cv_results$quantile_table_CV) <- c(paste0("X.", 1:length(X.dim)), "Accuracy", "Balanced Accuracy", "Precision", "Recall", "F1 score")
  return(cv_results)
}
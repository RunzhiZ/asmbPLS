#' Cross-validation for asmbPLS-DA to find the best combinations of quantiles for prediction
#'
#' Function to find the best combinations of quantiles used for prediction via
#' cross-validation. Usually should be conducted before 
#' \code{\link[asmbPLS]{asmbPLSDA.fit}} to obtain the quantile combinations.
#' 
#' @param X.matrix Predictors matrix. Samples in rows, variables in columns
#' @param Y.matrix Outcome matrix. Samples in rows, this is a matrix with one 
#' column (binary) or multiple columns (more than 2 levels, dummy variables).
#' @param PLS.comp Number of PLS components in asmbPLS-DA.
#' @param X.dim A vector containing the number of predictors in each block 
#' (ordered).
#' @param quantile.comb.table A matrix containing user-defined quantile 
#' combinations used for CV, whose column number equals to the 
#' number of blocks.
#' @param outcome.type The type of the outcome Y. "\code{binary}" for binary 
#' outcome, and "\code{morethan2levels}" for categorical outcome with more than 2 
#' levels.
#' @param Method Decision rule used for CV. For binary outcome, the 
#' methods include "\code{fixed_cutoff}", "\code{Euclidean_distance_X}" and
#' "\code{Mahalanobis_distance_X}". For categorical outcome with more than 2 
#' levels, the methods include "\code{Max_Y}", "\code{Euclidean_distance_X}",
#' "\code{Mahalanobis_distance_X}", "\code{Euclidean_distance_Y}", and 
#' "\code{PCA_Mahalanobis_distance_Y}".
#' @param k The number of folds of CV procedure. The default is 5.
#' @param ncv The number of repetitions of CV. The default is 5.
#' @param center A logical value indicating whether weighted mean center should 
#' be implemented for \code{X.matrix} and \code{Y.matrix}. The default is TRUE.
#' @param scale  A logical value indicating whether scale should be 
#' implemented for \code{X.matrix}. The default is TRUE.
#' 
#' @return 
#' \code{asmbPLSDA.cv} returns a list containing the following components:
#' \item{quantile_table_CV}{A matrix containing the selected quantile 
#' combination and the corresponding accuracy of CV for each PLS component.}
#' \item{CV_results}{A list containing the details of the CV results for each PLS 
#' component.}
#' \item{CV_index}{A list containing the validation and training index for each 
#' cross validation fold.}
#' 
#' @examples
#' ## Use the example dataset
#' data(asmbPLSDA.cv.example)
#' 
#' ## cv to find the best quantile combinations for model fitting (binary outcome)
#' cv.binary.results <- asmbPLSDA.cv(
#' X.matrix = asmbPLSDA.cv.example$X.matrix, 
#' Y.matrix = asmbPLSDA.cv.example$Y.matrix.binary, 
#' PLS.comp = asmbPLSDA.cv.example$PLS.comp, 
#' X.dim = asmbPLSDA.cv.example$X.dim, 
#' quantile.comb.table = asmbPLSDA.cv.example$quantile.comb.table, 
#' outcome.type = "binary", Method = "fixed_cutoff", k = 10, ncv = 5,
#' center = TRUE, scale = TRUE)
#' quantile.comb.binary <- cv.binary.results$quantile_table_CV[,1:2]
#' 
#' ## cv to find the best quantile combinations for model fitting 
#' ## (categorical outcome with more than 2 levels)
#' cv.morethan2levels.results <- asmbPLSDA.cv(
#' X.matrix = asmbPLSDA.cv.example$X.matrix, 
#' Y.matrix = asmbPLSDA.cv.example$Y.matrix.morethan2levels, 
#' PLS.comp = asmbPLSDA.cv.example$PLS.comp, 
#' X.dim = asmbPLSDA.cv.example$X.dim, 
#' quantile.comb.table = asmbPLSDA.cv.example$quantile.comb.table, 
#' outcome.type = "morethan2levels", Method = "Max_Y", k = 10, ncv = 5,
#' center = TRUE, scale = TRUE)
#' quantile.comb.morethan2levels <- cv.morethan2levels.results$quantile_table_CV[,1:2]
#'  
#' ## asmbPLSDA fit (binary outcome)
#' asmbPLSDA.binary.results <- asmbPLSDA.fit(
#' X.matrix = asmbPLSDA.cv.example$X.matrix, 
#' Y.matrix = asmbPLSDA.cv.example$Y.matrix.binary, 
#' PLS.comp = asmbPLSDA.cv.example$PLS.comp, 
#' X.dim = asmbPLSDA.cv.example$X.dim, 
#' quantile.comb = quantile.comb.binary,
#' "binary")
#' 
#' ## asmbPLSDA fit (categorical outcome with more than 2 levels)
#' asmbPLSDA.morethan2levels.results <- asmbPLSDA.fit(
#' X.matrix = asmbPLSDA.cv.example$X.matrix, 
#' Y.matrix = asmbPLSDA.cv.example$Y.matrix.morethan2levels, 
#' PLS.comp = asmbPLSDA.cv.example$PLS.comp, 
#' X.dim = asmbPLSDA.cv.example$X.dim, 
#' quantile.comb = quantile.comb.morethan2levels,
#' "morethan2levels")
#' 
#' @export
#' @useDynLib asmbPLS, .registration=TRUE
#' @importFrom Rcpp sourceCpp

asmbPLSDA.cv <- function(X.matrix, Y.matrix, PLS.comp, X.dim, quantile.comb.table, outcome.type, Method, k = 5, ncv = 5, center = TRUE, scale = TRUE) {
  ## error check
  stopifnot(!missing(X.matrix),
            !missing(Y.matrix),
            !missing(PLS.comp),
            !missing(X.dim),
            !missing(quantile.comb.table),
            !missing(outcome.type),
            !missing(Method),
            is.matrix(X.matrix),
            is.matrix(Y.matrix),
            is.matrix(quantile.comb.table),
            is.numeric(PLS.comp))
  return(asmbPLSDA_CV(X.matrix, Y.matrix, PLS.comp, X.dim, quantile.comb.table, outcome.type, Method, k, ncv, center, scale))
}
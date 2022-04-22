#' Cross-validation for asmbPLS to find the best combinations of quantiles for prediction
#'
#' Function to find the best combinations of quantiles used for prediction via
#' cross-validation. Usually should be conducted before 
#' \code{\link[asmbPLS]{asmbPLS.fit}} to obtain the quantile combinations.
#' 
#' @param X.matrix Predictors matrix. Samples in rows, variables in columns
#' @param Y.matrix Outcome matrix. Samples in rows, this is a matrix with one 
#' column (continuous variable). The outcome could be imputed survival time. 
#' For survival time with right-censored survival time and event indicator, the 
#' right censored time could be imputed by \code{\link{meanimp}}.
#' @param PLS.comp Number of PLS components in asmbPLS.
#' @param X.dim A vector containing the number of predictors in each block 
#' (ordered).
#' @param quantile.comb.table A matrix containing user-defined quantile 
#' combinations used for CV, whose column number equals to the 
#' number of blocks.
#' @param k The number of folds of CV procedure. The default is 5.
#' @param Y.indicator A vector containing the event indicator for each sample, 
#' whose length is equal to the number of samples. This vector allows the ratio 
#' of observed/unobserved to be the same in the training set and validation set. 
#' Observed = 1, and unobserved = 0. If other types of outcome data rather than 
#' survival outcome is used, you can use a vector with all components = 1 
#' instead.
#' @param only.observe Whether only observed samples in the validation set 
#' should be used for calculating the MSE for CV. The default is TRUE.
#' 
#' @return 
#' \code{asmbPLS.cv} returns a list containing the following components:
#' \item{quantile.comb}{A matrix containing the selected quantile 
#' combination and the corresponding MSE of CV for each PLS component.}
#' \item{time}{Computation time.}
#' 
#' @examples
#' ## Use the example dataset
#' data(asmbPLS.cv.example)
#' 
#' ## cv to find the best quantile combinations for model fitting
#' cv.results <- asmbPLS.cv(X.matrix = X.matrix, Y.matrix = Y.matrix, PLS.comp 
#' = PLS.comp, X.dim = X.dim, quantile.comb.table = quantile.comb.table, k = 5, 
#' Y.indicator = Y.indicator, only.observe = T)
#' quantile.comb <- cv.results$quantile.comb[,1:2]
#'  
#' ## asmbPLS fit
#' asmbPLS.results <- asmbPLS.fit(X.matrix = X.matrix, Y.matrix = Y.matrix, 
#' PLS.comp = PLS.comp, X.dim = X.dim, quantile.comb = quantile.comb)
#' 
#' @export
#' @useDynLib asmbPLS, .registration=TRUE
#' @importFrom Rcpp sourceCpp

asmbPLS.cv <- function(X.matrix, Y.matrix, PLS.comp, X.dim, quantile.comb.table, k = 5, Y.indicator, only.observe = T) {
  ## error check
  stopifnot(!missing(X.matrix),
            !missing(Y.matrix),
            !missing(PLS.comp),
            !missing(X.dim),
            !missing(quantile.comb.table),
            !missing(Y.indicator),
            is.matrix(X.matrix),
            is.matrix(Y.matrix),
            is.matrix(quantile.comb.table))
  ## cv procedure
  time_start <- Sys.time()
  PLS_table <- matrix(NA, nrow = PLS.comp, ncol = length(X.dim))
  quantile.comb <- matrix(NA, nrow = PLS.comp, ncol = length(X.dim)+1)
  n_CV <- nrow(X.matrix)
  all_obs_indicator <- (!any(Y.indicator== 0))
  if (all_obs_indicator) {
    Y_observed_index <- cbind(which(Y.indicator == 1), cut(1:table(Y.indicator)["1"], k, labels = F))
  } else {
    Y_censored_index <- cbind(which(Y.indicator == 0), cut(1:table(Y.indicator)["0"], k, labels = F))
    Y_observed_index <- cbind(which(Y.indicator == 1), cut(1:table(Y.indicator)["1"], k, labels = F))
  }
  cv.results <- matrix(NA, nrow = nrow(quantile.comb.table), ncol = k)
  for (PLS_term in 1:PLS.comp) {
    for (g in 1:k) {
      if (all_obs_indicator) {
        validation_index <- Y_observed_index[which(Y_observed_index[, 2] == g), 1]
        training_index <- Y_observed_index[which(Y_observed_index[, 2] != g), 1]
      } else {
        if (only.observe) {
          validation_index <- Y_observed_index[which(Y_observed_index[, 2] == g), 1]
          training_index <- c(Y_censored_index[which(Y_censored_index[, 2] != g), 1], Y_observed_index[which(Y_observed_index[, 2] != g), 1])
        } else {
          validation_index <- c(Y_censored_index[which(Y_censored_index[, 2] == g), 1], Y_observed_index[which(Y_observed_index[, 2] == g), 1])
          training_index <- c(Y_censored_index[which(Y_censored_index[, 2] != g), 1], Y_observed_index[which(Y_observed_index[, 2] != g), 1])
        }
      }
      X_train_CV <- X.matrix[training_index, ]
      X_validation_CV <- X.matrix[validation_index, ]
      Y_train_CV <- as.matrix((Y.matrix[training_index, ]))
      Y_validation_CV <- as.matrix((Y.matrix[validation_index, ]))
      for (m in 1:nrow(quantile.comb.table)) {
        PLS_table[PLS_term, ] <- quantile.comb.table[m, ]
        asmbPLS_results_CV <- asmbPLS_fit_rcpp(X_train_CV, Y_train_CV, PLS_term, X.dim, matrix(PLS_table[1:PLS_term, ],ncol = length(X.dim)))
        Y_asmbPLS_CV <- asmbPLS_predict_rcpp(X_validation_CV, PLS_term, asmbPLS_results_CV)
        cv.results[m, g] <- result_compare(Y_asmbPLS_CV, Y_validation_CV)[[1]]
      }
    }
    cv.results.mean <- apply(cv.results,1,mean)
    PLS_table[PLS_term, ] <- quantile.comb[PLS_term, 1:length(X.dim)] <- quantile.comb.table[which.min(cv.results.mean), ]
    quantile.comb[PLS_term, length(X.dim)+1] <- cv.results.mean[which.min(cv.results.mean)]
  }
  time_end <- Sys.time()
  time_diff <- difftime(time_end, time_start, units = "secs")
  print(paste0("Time for asmbPLS CV: ",round(time_diff,2)))
  return(list(quantile.comb = quantile.comb, 
              time = time_diff))
}
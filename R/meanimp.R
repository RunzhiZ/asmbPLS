#' Mean imputation for the survival time
#'
#' In this approach, \eqn{\mu} can be computed using the familiar sample mean
#' formula provided the censored values are imputed.
#' 
#' @param survival.data A matrix of two columns with the first column indicates
#' the survival time and the second column indicates the event indicator (1 and 
#' 0, where 1 indicates observed event and 0 indicates unobserved event).
#' @param round Whether survival time should be rounded, default = FALSE.
#' 
#' @return 
#' \code{meanimp} returns a list containing the following components:
#' \item{imputed_table}{A matrix containing the original survival data and the 
#' imputed time.}
#' \item{KM_table}{Kaplan-Meier estimator of failure times.}
#' 
#' @examples
#' ## Generate the survival data
#' data_test <- matrix(c(1, 1, 1, 2.5, 5, 7, 1, 1, 0, 1, 0, 1), ncol = 2)
#' 
#' ## Mean imputation
#' meanimp(data_test, round = FALSE)
#' 
#' @export
#' @useDynLib asmbPLS, .registration=TRUE
#' @importFrom Rcpp sourceCpp

meanimp <- function(survival.data, round = FALSE) {
  colnames(survival.data) <- c("Survival_time", "Event_indicator")
  n <- nrow(survival.data)
  survival.data <- cbind(survival.data, 1:n, NA)
  colnames(survival.data)[3:4] <- c("ID", "Imputed_time")
  survival.data <- survival.data[order(survival.data[, "Survival_time"]), ]
  survival.data[n, "Event_indicator"]<-1
  if (round == T) {
    survival.data[, "Survival_time"] <- round(survival.data[, "Survival_time"])
  }
  KM_table <- KM.estimator(survival.data)
  
  censored_index <- which(survival.data[, "Event_indicator"] == 0)
  
  for (i in 1:length(censored_index)) {
    censored_time <- survival.data[censored_index[i], "Survival_time"]
    min_survival_time_index <- min(which(KM_table[, "Time"] > censored_time))
    numerator_table <- matrix(NA, nrow = nrow(KM_table) - min_survival_time_index + 1, ncol = 2)
    colnames(numerator_table) <- c("Tao", "Survival_diff")
    numerator_table[, "Tao"] <- KM_table[min_survival_time_index:nrow(KM_table), "Time"]
    if (min_survival_time_index == 1) {
      numerator_table[, "Survival_diff"] <- c(1, KM_table[min_survival_time_index:(nrow(KM_table)-1), "Survival_total"]) -
        KM_table[min_survival_time_index:nrow(KM_table), "Survival_total"]
    } else {
      numerator_table[, "Survival_diff"] <- KM_table[(min_survival_time_index-1):(nrow(KM_table)-1), "Survival_total"] -
        KM_table[min_survival_time_index:nrow(KM_table), "Survival_total"]
    }
    numerator <- sum(numerator_table[, "Tao"] * numerator_table[, "Survival_diff"])
    if (censored_time >= min(KM_table[, "Time"])) {
      denominator <- KM_table[max(which(KM_table[, "Time"] <= censored_time)), "Survival_total"]
    } else {
      denominator <- 1
    }
    survival.data[censored_index[i], "Imputed_time"] <- numerator / denominator
  }
  
  survival.data <- survival.data[order(survival.data[, "ID"]), ]
  survival.data[which(survival.data[,"Event_indicator"] == 1),"Imputed_time"] <- survival.data[which(survival.data[,"Event_indicator"] == 1),"Survival_time"]
  return(list(imputed_table = survival.data,
              KM_table = KM_table))
}
#' Create the quantile combination set for asmbPLS
#'
#' Create the quantile combination set given quantile set for each block
#' 
#' @param quantile.list A list containing the quantile set for each block.
#' 
#' @return 
#' \item{quantile.comb}{The quantile combination used for asmbPLS}
#' 
#' @examples
#' ## Generate quantile set for each block
#' ## For example, we have three blocks
#' quantile_1 <- c(0.999, 0.9992, 0.9994, 0.9996, 0.9998)
#' quantile_2 <- c(0.96, 0.97, 0.98, 0.99, 0.995)
#' quantile_3 <- c(0.95, 0.96, 0.97, 0.98, 0.99)
#' quantilelist <- list(quantile_1, quantile_2, quantile_3)
#' quantile.comb <- quantile.combination(quantilelist)
#' 
#' @export
#' @useDynLib asmbPLS, .registration=TRUE
#' @importFrom Rcpp sourceCpp

quantile.combination <- function(quantilelist) {
  n <- length(quantilelist)
  eval(parse(text = paste0("output <- expand.grid(", paste(paste0("quantilelist[[", 1:n, "]]"), collapse = ", "), ")")))
  colnames(output) <- paste0("block_", 1:n)
  output <- as.matrix(output)
  return(output)
}
#' Converts a class vector to binary class matrix
#'
#' Converts a class vector to binary class matrix with the number of columns equal to the number of levels  
#' 
#' @param categorical.vector A class vector
#' 
#' @return 
#' \item{output.matrix}{The output binary class matrix, which can be used for asmbPLS-DA.}
#' 
#' @examples
#' ## Generate a class vector
#' vector.test <- factor(c(1,1,1,2,2,2,3,3,3), levels = c(1,2,3))
#' 
#' ## convert
#' output.matrix <- to.categotical(vector.test)
#' 
#' @export
#' @useDynLib asmbPLS, .registration=TRUE
#' @importFrom Rcpp sourceCpp

to.categorical <- function(data) {
  if(is.factor(data)) {
    data_u <- levels(data)
  } else {
    data_u <- unique(data)
  }
  output_matrix <- matrix(0, nrow = length(data), ncol = length(data_u))
  colnames(output_matrix) <- data_u
  for(i in 1:length(data_u)) {
    output_matrix[which(data == data_u[i]), i] <- 1
  }
  return(output_matrix)
}
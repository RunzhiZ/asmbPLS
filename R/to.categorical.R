#' Converts a class vector to a binary class matrix
#'
#' This function converts a class vector to a binary class matrix, with the number of columns equal to the number of levels in the input vector. Each row of the output matrix corresponds to a single observation in the input vector, and the columns represent the different classes in the input vector. A value of 1 in a particular column indicates that the corresponding observation belongs to that class, while a value of 0 indicates that it does not.
#'
#' @param categorical.vector A factor or character vector representing the class labels.
#' 
#' @return A binary class matrix with the number of rows equal to the length of the input vector, and the number of columns equal to the number of unique levels in the input vector. The row and column names of the output matrix are set to the levels of the input vector.
#' 
#' @examples
#' ## Generate a class vector
#' vector.test <- factor(c(1,1,1,2,2,2,3,3,3), levels = c(1,2,3))
#' 
#' ## Convert the class vector to binary class matrix
#' output.matrix <- to.categorical(vector.test)
#' 
#' @export
#' @useDynLib asmbPLS, .registration=TRUE

to.categorical <- function(categorical.vector) {
  if(is.factor(categorical.vector)) {
    data_u <- levels(categorical.vector)
  } else {
    data_u <- unique(categorical.vector)
  }
  output_matrix <- matrix(0, nrow = length(categorical.vector), ncol = length(data_u))
  colnames(output_matrix) <- data_u
  for(i in 1:length(data_u)) {
    output_matrix[which(categorical.vector == data_u[i]), i] <- 1
  }
  return(output_matrix)
}
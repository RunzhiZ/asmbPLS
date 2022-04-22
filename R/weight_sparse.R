weight_sparse <- function(x, lambda){
  if(lambda == max(abs(x))) {
    uniq_x <- unique(unlist(x))
    lambda <- sort(uniq_x)[length(uniq_x)-1]
  }
  sparse_w = matrix(apply(x, 1, function(x, lambda) {return(sign(x)*max(abs(x)-lambda, 0))}, lambda), nrow = nrow(x))
  return(sparse_w)
}

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
arma::vec quantile_f(arma::vec V, 
                     arma::vec P) {
  return quantile(V, P);
}
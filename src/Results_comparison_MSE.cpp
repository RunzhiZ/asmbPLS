#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
double Results_comparison_MSE(arma::mat Y_predict, 
                              arma::mat Y_true) {
  double n_sample = Y_true.n_rows;
  arma::mat Y_diff = Y_predict - Y_true;
  double MSE = as_scalar(Y_diff.t()*Y_diff/n_sample);
  return(MSE);
} 
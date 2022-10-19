#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
double Results_comparison_accuracy(arma::mat Y_predict, 
                                   arma::mat Y_true) {
  int Y_col = Y_true.n_cols;
  double n_sample = Y_true.n_rows;
  double n_match = 0;
  if (Y_col == 1) {
    for (int i = 0; i < n_sample; ++i) {
      if (Y_predict(i, 0) == Y_true(i, 0)) {
        n_match = n_match + 1;
      }
    }
  } else {
    for (int i = 0; i < n_sample; ++i) {
      bool check_equal = approx_equal(Y_predict.row(i), Y_true.row(i), "absdiff", 0.00001);
      if (check_equal) {
        n_match = n_match + 1;
      }
    }
  }
  double accuracy = n_match/n_sample;
  return(accuracy);
} 
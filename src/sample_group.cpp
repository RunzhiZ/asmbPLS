#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector sample_group(double n, 
                           double K_input) {
  NumericVector output(n);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < K_input; ++j) {
      double temp = n/K_input;
      if((i >= j * temp) & (i < (j + 1) * temp + 1)) {
        output(i) = j;
      }
    }
  }
  return(output);
}
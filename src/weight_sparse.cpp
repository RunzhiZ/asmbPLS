#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
arma::colvec weight_sparse(arma::colvec input, 
                           double lambda) {
  if(lambda == max(abs(input))){
    arma::colvec x_unique = unique(input);
    arma::colvec x_sort = sort(x_unique, "descend");
    lambda = x_sort.at(1);
  }
  arma::colvec part_1 = abs(input) - lambda;
  int n_part_1 = part_1.size();
  arma::colvec part_max(n_part_1);
  for (int i = 0; i < n_part_1; ++i){
    part_max.at(i) = std::max(part_1.at(i), 0.0);
  }
  arma::colvec part_sign = sign(input);
  arma::colvec final_output(n_part_1);
  for (int i = 0; i < n_part_1; ++i){
    final_output.at(i) = part_max.at(i)*part_sign.at(i);
  }
  return(final_output);
}
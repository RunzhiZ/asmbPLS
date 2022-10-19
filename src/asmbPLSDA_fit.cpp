#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
List asmbPLSDA_fit(arma::mat X_matrix, 
                   arma::mat Y_matrix, 
                   int PLS_term, 
                   NumericVector X_dim, 
                   arma::mat percent, 
                   String outcome_type,
                   Nullable<LogicalVector> center = R_NilValue, 
                   Nullable<LogicalVector> scale = R_NilValue) {
  
  Function asmbPLSDA_binary_fit = Environment::namespace_env("asmbPLS")["asmbPLSDA_binary_fit"];
  Function asmbPLSDA_morethantwo_fit = Environment::namespace_env("asmbPLS")["asmbPLSDA_morethantwo_fit"];
  // Function asmbPLSDA_binary_fit("asmbPLSDA_binary_fit");
  // Function asmbPLSDA_morethantwo_fit("asmbPLSDA_morethantwo_fit");
  
  List output;
  if (outcome_type == "binary") {
    output = asmbPLSDA_binary_fit(X_matrix, Y_matrix, PLS_term, X_dim, percent, center, scale);
  }
  if (outcome_type == "morethan2levels") {
    output = asmbPLSDA_morethantwo_fit(X_matrix, Y_matrix, PLS_term, X_dim, percent, center, scale);
  }
  return(output);
}




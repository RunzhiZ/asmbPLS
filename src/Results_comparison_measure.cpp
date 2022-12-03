#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
arma::rowvec Results_comparison_measure(arma::mat Y_predict, 
                                        arma::mat Y_true) {
  double Y_col = Y_true.n_cols;
  double n_sample = Y_true.n_rows;
  double n_match;
  double n_TP;
  double n_FP;
  double n_FN;
  double accuracy;
  double precision;
  double recall;
  double F1;
  arma::mat output_temp(Y_col, 4);
  
  for (int i = 0; i < Y_col; ++i) {
    arma::uvec temp_accu = find(Y_predict.col(i) == Y_true.col(i));
    arma::uvec temp_TP = find(Y_predict.col(i) == 1 && Y_true.col(i) == 1);
    arma::uvec temp_TN = find(Y_predict.col(i) == 0 && Y_true.col(i) == 0);
    arma::uvec temp_FP = find(Y_predict.col(i) == 1 && Y_true.col(i) == 0);
    arma::uvec temp_FN = find(Y_predict.col(i) == 0 && Y_true.col(i) == 1);
    
    n_match = temp_accu.n_rows;
    n_TP = temp_TP.n_rows;
    n_FP = temp_FP.n_rows;
    n_FN = temp_FN.n_rows;
    
    accuracy = n_match/n_sample;
    precision = n_TP/(n_TP + n_FN);
    recall = n_TP/(n_TP + n_FP);
    F1 = 2 * precision * recall/(precision + recall);
    
    output_temp(i, 0) = accuracy;
    output_temp(i, 1) = precision;
    output_temp(i, 2) = recall;
    output_temp(i, 3) = F1;
  }
  
  arma::rowvec output = mean(output_temp, 0);
  
  return(output);
}
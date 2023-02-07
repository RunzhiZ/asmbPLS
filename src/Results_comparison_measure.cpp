#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
arma::rowvec Results_comparison_measure(arma::mat Y_predict, 
                                        arma::mat Y_true, 
                                        String outcome_type) {
  double Y_col = Y_true.n_cols; // Number of groups of samples
  double n_match = 0;
  double n_TP = 0; // Number of true positive
  double n_TN = 0; // Number of true negative
  double n_FP = 0; // Number of false positive
  double n_FN = 0; // Number of false negative
  double balanced_accuracy_multicalss = 0;
  arma::rowvec output(1, 5);
  double n_recall_multiclass;
  
  for (int i = 0; i < Y_col; ++i) {
    arma::uvec temp_accu = find(Y_predict.col(i) == Y_true.col(i));
    arma::uvec temp_TP = find(Y_predict.col(i) == 1 && Y_true.col(i) == 1);
    arma::uvec temp_TN = find(Y_predict.col(i) == 0 && Y_true.col(i) == 0);
    arma::uvec temp_FP = find(Y_predict.col(i) == 1 && Y_true.col(i) == 0);
    arma::uvec temp_FN = find(Y_predict.col(i) == 0 && Y_true.col(i) == 1);
    
    double n_temp_TP = temp_TP.n_rows;
    double n_temp_TN = temp_TN.n_rows;
    double n_temp_FP = temp_FP.n_rows;
    double n_temp_FN = temp_FN.n_rows;
    
    n_match = n_match + temp_accu.n_rows;
    n_TP = n_TP + n_temp_TP;
    n_TN = n_TN + n_temp_TN;
    n_FP = n_FP + n_temp_FP;
    n_FN = n_FN + n_temp_FN;
    
    n_recall_multiclass = n_temp_TP/(n_temp_TP + n_temp_FN);
    balanced_accuracy_multicalss = balanced_accuracy_multicalss + n_recall_multiclass;
  }
  double accuracy = (n_TP + n_TN)/(n_TP + n_TN + n_FP + n_FN);
  double precision = n_TP/(n_TP + n_FP);
  double recall = n_TP/(n_TP + n_FN); // Sensitivity also
  double specificity = n_TN/(n_TN + n_FP);
  double F1 = 2 * precision * recall/(precision + recall);
  
  double balanced_accuracy = 0;
  double balanced_accuracy_binary = (recall + specificity)/2;
  balanced_accuracy_multicalss = balanced_accuracy_multicalss/Y_col;
  if (outcome_type == "binary") {
    balanced_accuracy = balanced_accuracy_binary;
  }
  if (outcome_type == "multiclass") {
    balanced_accuracy = balanced_accuracy_multicalss;
  }
  
  output(0) = accuracy;
  output(1) = balanced_accuracy;
  output(2) = precision;
  output(3) = recall;
  output(4) = F1;
  
  return(output);
}
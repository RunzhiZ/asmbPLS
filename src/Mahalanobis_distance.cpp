#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
arma::mat Mahalanobis_distance(arma::mat matrix_fit, 
                               arma::mat matrix_predict, 
                               arma::mat F_matrix, 
                               String outcome_type) {
  
  double n_row_predict = matrix_predict.n_rows; // Number of samples in new data
  double n_col_predict = matrix_predict.n_cols; // Number of features
  double n_col_F = F_matrix.n_cols; // Number of groups of samples
  arma::mat predict_output(n_row_predict, n_col_F, arma::fill::zeros); // Classification output
  
  if (outcome_type == "binary") {
    arma::mat matrix_fit_1 = matrix_fit.rows(find(F_matrix == 0));
    arma::mat matrix_fit_2 = matrix_fit.rows(find(F_matrix == 1));
    arma::rowvec g1_centroid = mean(matrix_fit_1, 0);
    arma::rowvec g2_centroid = mean(matrix_fit_2, 0);
    // S_inverse matrix
    arma::mat g1_S = cov(matrix_fit_1);
    arma::mat g2_S = cov(matrix_fit_2);
    double n_1 = matrix_fit_1.n_rows;
    double n_2 = matrix_fit_2.n_rows;
    arma::mat S = ((n_1 - 1)*g1_S + (n_2 - 1)*g2_S)/(n_1 + n_2 -2); 
    arma::mat S_inv = inv(S);
    
    arma::mat matrix_predict_1 = matrix_predict.each_row() - g1_centroid;
    arma::mat matrix_predict_2 = matrix_predict.each_row() - g2_centroid;
    arma::mat matrix_predict_1_temp = sqrt(matrix_predict_1*S_inv*matrix_predict_1.t());
    arma::mat matrix_predict_2_temp = sqrt(matrix_predict_2*S_inv*matrix_predict_2.t());
    predict_output.elem(find(matrix_predict_2_temp.diag() < matrix_predict_1_temp.diag())).fill(1);
  }

  if (outcome_type == "multiclass") {
    arma::mat g_score(n_row_predict, n_col_F);
    arma::mat S(n_col_predict, n_col_predict, arma::fill::zeros); // S matrix
    double n = 0;
    for (int i = 0; i < n_col_F; ++i) {
      arma::mat matrix_fit_temp = matrix_fit.rows(find(F_matrix.col(i) == 1));
      arma::mat S_g = cov(matrix_fit_temp);
      double n_g = matrix_fit_temp.n_rows;
      n = n + n_g;
      S = S + (n_g - 1)*S_g;
    }
    // S_inverse matrix
    S = S/(n - n_col_F);
    arma::mat S_inv = inv(S); // S_inverse matrix
    
    for (int i = 0; i < n_col_F; ++i) {
      arma::mat matrix_fit_temp = matrix_fit.rows(find(F_matrix.col(i) == 1));
      arma::rowvec g_centroid_temp = mean(matrix_fit_temp, 0);
      arma::mat matrix_predict_g = matrix_predict.each_row() - g_centroid_temp;
      arma::mat matrix_predict_temp = sqrt(matrix_predict_g*S_inv*matrix_predict_g.t());
      g_score.col(i) = matrix_predict_temp.diag();
    }
    arma::ucolvec minvec = arma::index_min(g_score, 1);
    for (int i = 0; i < n_row_predict; ++i) {
      predict_output(i, minvec(i)) = 1;
    }
  }
  return(predict_output);
} 
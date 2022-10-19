#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
arma::mat Euclidean_distance(arma::mat matrix_fit, 
                             arma::mat matrix_predict, 
                             arma::mat F_matrix, 
                             String outcome_type) {
  
  int n_row_predict = matrix_predict.n_rows;
  int n_col_F = F_matrix.n_cols;
  arma::mat predict_output(n_row_predict, n_col_F);
  
  if (outcome_type == "binary") {
    arma::colvec F_group = unique(F_matrix);
    arma::uvec g1_index = find(F_matrix == F_group[0]);
    arma::uvec g2_index = find(F_matrix == F_group[1]);
    arma::mat matrix_fit_1 = matrix_fit.rows(g1_index);
    arma::mat matrix_fit_2 = matrix_fit.rows(g2_index);
    arma::colvec g1_centroid = mean(matrix_fit_1, 0).t();
    arma::colvec g2_centroid = mean(matrix_fit_2, 0).t();
    for (int i = 0; i < n_row_predict; ++i) {
      arma::colvec predict_temp = matrix_predict.row(i).t();
      arma::colvec predict_temp_1 = predict_temp - g1_centroid;
      arma::colvec predict_temp_2 = predict_temp - g2_centroid;
      double g1_score = sqrt(as_scalar(predict_temp_1.t()*predict_temp_1));
      double g2_score = sqrt(as_scalar(predict_temp_2.t()*predict_temp_2));
      if (g1_score < g2_score) {
        predict_output(i, 0) = F_group[0];
      } else {
        predict_output(i, 0) = F_group[1];
      }
    }
  }
  
  // matrix_fit_group
  List matrix_fit_group(n_col_F);
  List g_centroid(n_col_F);
  if (outcome_type == "morethan2levels") {
    predict_output = predict_output.zeros();
    for (int i = 0; i < n_col_F; ++i) {
      arma::uvec temp = find(F_matrix.col(i) == 1);
      arma::mat matrix_fit_temp = matrix_fit.rows(temp);
      matrix_fit_group[i] = matrix_fit_temp;
      arma::colvec g_centroid_temp = mean(matrix_fit_temp, 0).t();
      g_centroid[i] = g_centroid_temp;
    }
    for (int i = 0; i < n_row_predict; ++i) {
      arma::colvec g_score(n_col_F, arma::fill::zeros);
      for (int j = 0; j < n_col_F; ++j) {
        arma::colvec g_centroid_temp = g_centroid[j];
        arma::colvec predict_temp = matrix_predict.row(i).t();
        arma::colvec predict_temp_g = predict_temp - g_centroid_temp;
        double g_score_temp = sqrt(as_scalar(predict_temp_g.t()*predict_temp_g));
        g_score(j) = g_score_temp;
      } 
      int index_min = g_score.index_min();
      predict_output(i, index_min) = 1;
    }
  }
  
  return(predict_output);
}
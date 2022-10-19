#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
arma::mat PCA_Mahalanobis_distance(arma::mat Y_fit, 
                                   arma::mat Y_predict) {
  
  // call PCA function
  Environment stats("package:stats");
  Function prcomp = stats["prcomp"];
  
  int Y_col = Y_fit.n_cols;
  int n_predict = Y_predict.n_rows;
  arma::mat E(Y_col, Y_col, arma::fill::eye);
  List PCA_output = prcomp(Y_fit);
  arma::mat PCA_T = PCA_output["x"];
  arma::mat PCA_P = PCA_output["rotation"]                             ;
  arma::rowvec PCA_mean = PCA_output["center"];
  arma::mat mean_matrix(n_predict, Y_col);
  arma::mat mean_matrix_C(Y_col, Y_col);
  arma::mat PCA_T_1 = PCA_T.cols(0, Y_col - 2);
  arma::mat PCA_P_1 = PCA_P.cols(0, Y_col - 2);
  arma::mat S_T = PCA_T_1.t()*PCA_T_1;
  arma::mat S_T_inv = inv(S_T);
  
  for (int i = 0; i < Y_col; ++i) {
    for (int j = 0; j < n_predict; ++j) {
      mean_matrix(j, i) = arma::as_scalar(PCA_mean(i));
    }
    for (int j = 0; j < Y_col; ++j) {
      mean_matrix_C(j, i) = arma::as_scalar(PCA_mean(i));
    }
  }
  
  arma::mat T_predict = (Y_predict - mean_matrix)*PCA_P_1;
  arma::mat C = (E - mean_matrix_C)*PCA_P_1;
  arma::mat predict_output(n_predict, Y_col);
  
  for (int i = 0; i < n_predict; ++i) {
    arma::colvec g_score(Y_col, arma::fill::zeros);
    for (int j = 0; j < Y_col; ++j) {
      arma::colvec g_centroid_temp = C.row(j).t();
      arma::colvec predict_temp = T_predict.row(i).t();
      arma::colvec predict_temp_g = predict_temp - g_centroid_temp;
      double g_score_temp = sqrt(as_scalar(predict_temp_g.t()*S_T_inv*predict_temp_g));
      g_score(j) = g_score_temp;
    } 
    int index_min = g_score.index_min();
    predict_output(i, index_min) = 1;
  }
  
  return(predict_output);
}
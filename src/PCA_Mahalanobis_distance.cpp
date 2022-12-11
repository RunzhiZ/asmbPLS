#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
arma::mat PCA_Mahalanobis_distance(arma::mat Y_fit, 
                                   arma::mat Y_predict) {
  
  // call PCA function
  Environment stats("package:stats");
  Function prcomp = stats["prcomp"];
  
  int Y_col = Y_fit.n_cols; // Number of groups of features 
  int n_row_predict = Y_predict.n_rows; // Number of samples in new data
  arma::mat predict_output(n_row_predict, Y_col); // Classification output
  
  arma::mat E(Y_col, Y_col, arma::fill::eye);
  List PCA_output = prcomp(Y_fit); // Results of PCA
  arma::mat PCA_T = PCA_output["x"]; // PCA score
  arma::mat PCA_P = PCA_output["rotation"]; // PCA loading                             ;
  arma::rowvec PCA_mean = PCA_output["center"]; // Column mean of Y
  
  arma::mat mean_matrix(n_row_predict, Y_col);
  arma::mat mean_matrix_C(Y_col, Y_col);
  arma::mat PCA_T_1 = PCA_T.cols(0, Y_col - 2);
  arma::mat PCA_P_1 = PCA_P.cols(0, Y_col - 2);
  arma::mat S_T = PCA_T_1.t()*PCA_T_1;
  arma::mat S_T_inv = inv(S_T);
  
  mean_matrix.each_row() -= PCA_mean;
  mean_matrix_C.each_row() -= PCA_mean;
  
  arma::mat T_predict = (Y_predict - mean_matrix)*PCA_P_1;
  arma::mat C = (E - mean_matrix_C)*PCA_P_1;
  
  arma::mat g_score(n_row_predict, Y_col);
  
  for (int i = 0; i < Y_col; ++i) {
    arma::rowvec g_centroid_temp = C.row(i);
    arma::mat T_predict_g = T_predict.each_row() - g_centroid_temp;
    arma::mat T_predict_temp = sqrt(T_predict_g*S_T_inv*T_predict_g.t());
    g_score.col(i) = T_predict_temp.diag();
  }
  arma::ucolvec minvec = arma::index_min(g_score, 1);
  for (int i = 0; i < n_row_predict; ++i) {
    predict_output(i, minvec(i)) = 1;
  }
  return(predict_output);
}
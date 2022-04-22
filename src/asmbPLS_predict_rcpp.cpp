#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
arma::mat asmbPLS_predict_rcpp(arma::mat newdata, int PLS_term_selected, List asmbPLS_results) {
  NumericVector X_dim = asmbPLS_results["X_dim"];
  
  // newdata scale
  int E_col = newdata.n_cols;
  arma::rowvec col_mean = as<arma::rowvec>(asmbPLS_results["col_mean"]);
  arma::rowvec col_sd = as<arma::rowvec>(asmbPLS_results["col_sd"]);
  for (int i = 0; i < E_col; ++i) {
    arma::mat temp = (newdata.cols(i, i) - arma::as_scalar(col_mean.col(i)))/arma::as_scalar(col_sd.col(i));
    newdata.cols(i, i) = temp;
  }
  double Y_mean = as<double>(asmbPLS_results["Y_mean"]);
  double Y_sd = as<double>(asmbPLS_results["Y_sd"]);
  
  // variables for convenient
  int B = X_dim.length();
  X_dim.push_front(0);
  int n_row = newdata.n_rows;
  arma::rowvec y_weight = as<arma::rowvec>(asmbPLS_results["y_weight"]);
  arma::mat Y_pred(n_row, PLS_term_selected);
  
  // temp variables
  arma::mat X_matrix_temp;
  arma::colvec w_temp;
  arma::colvec t_temp(n_row);
  arma::mat t_cbind(n_row, B);
  arma::colvec w_T_temp(B);
  arma::colvec t_T(n_row);
  arma::colvec p_temp;
  arma::colvec q;
  
  
  arma::mat x_score_mat(n_row*B, PLS_term_selected);
  
  for (int i = 0; i < PLS_term_selected; ++i) {
    for (int j = 0; j < B; ++j) {
      X_matrix_temp = newdata.cols(sum(X_dim[Range(0, j)]), sum(X_dim[Range(0, j+1)]) - 1);
      w_temp = as<arma::mat>(as<List>(asmbPLS_results["x_weight"])[j]).col(i);
      t_temp = X_matrix_temp*w_temp/(sqrt(X_dim[j+1]));
      x_score_mat.submat((j)*n_row, i, (j+1)*n_row-1, i) = t_temp;
      t_cbind.col(j) = t_temp;
    }
    w_T_temp = as<arma::mat>(asmbPLS_results["x_super_weight"]).col(i);
    t_T = t_cbind*w_T_temp;
    for (int j = 0; j < B; ++j) {
      X_matrix_temp = newdata.cols(sum(X_dim[Range(0, j)]), sum(X_dim[Range(0, j+1)]) - 1);
      p_temp = as<arma::mat>(as<List>(asmbPLS_results["x_loading"])[j]).col(i);
      newdata.cols(sum(X_dim[Range(0, j)]), sum(X_dim[Range(0, j+1)]) - 1) = X_matrix_temp - t_T*p_temp.t();
    }
    q = y_weight.col(i);
    Y_pred.col(i) = t_T*q.t();
  }
  arma::mat Y_pred_final = sum(Y_pred, 1);
  arma::mat temp = Y_pred_final*Y_sd + Y_mean;
  Y_pred_final = temp;
  return (Y_pred_final);
}


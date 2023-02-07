#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
List asmbPLS_predict(List asmbPLS_results,
                     arma::mat newdata, 
                     int PLS_term_selected) {
  
  NumericVector X_dim = asmbPLS_results["X_dim"];
  LogicalVector center = asmbPLS_results["center"];
  LogicalVector scale = asmbPLS_results["scale"];
  arma::rowvec X_col_mean = as<arma::rowvec>(asmbPLS_results["X_col_mean"]);
  arma::rowvec X_col_sd = as<arma::rowvec>(asmbPLS_results["X_col_sd"]);
  double Y_col_mean = as<double>(asmbPLS_results["Y_col_mean"]);
  double Y_col_sd = as<double>(asmbPLS_results["Y_col_sd"]);
  List x_weight = as<List>(asmbPLS_results["X_weight"]);
  arma::mat x_super_weight = as<arma::mat>(asmbPLS_results["X_super_weight"]);
  List x_loading = as<List>(asmbPLS_results["X_loading"]);
  arma::rowvec y_weight = as<arma::rowvec>(asmbPLS_results["Y_weight"]);
  
  int E_col = newdata.n_cols; // Number of features
  int B = X_dim.length(); // Number of blocks
  int F_col = 1; // Number of columns of Y.matrix
  int n_row = newdata.n_rows; // Number of samples for new data
  
  // Newdata scale
  if (center[0]) {
    if (scale[0]) {
      // center = 1 and scale = 1
      for (int i = 0; i < E_col; ++i) {
        arma::colvec temp = (newdata.col(i) - arma::as_scalar(X_col_mean(i)))/arma::as_scalar(X_col_sd(i));
        newdata.col(i) = temp;
      }
    } else {
      // center = 1 and scale = 0
      for (int i = 0; i < E_col; ++i) {
        arma::colvec temp = (newdata.col(i) - arma::as_scalar(X_col_mean(i)));
        newdata.col(i) = temp;
      }
    }
  } else {
    if (scale[0]) {
      // center = 0 and scale = 1
      for (int i = 0; i < E_col; ++i) {
        arma::colvec temp = newdata.col(i)/arma::as_scalar(X_col_sd(i));
        newdata.col(i) = temp;
      }
    }
  } // center = 0 and scale = 1, no action 
  
  X_dim.push_front(0);
  
  arma::mat Y_pred(n_row, F_col, arma::fill::zeros); // Prediction for new data
  arma::mat t_T_all(n_row, PLS_term_selected); // Super score for new data
  
  // Temp variables
  arma::mat X_matrix_temp;
  arma::colvec w_temp;
  arma::colvec t_temp(n_row);
  arma::mat t_cbind(n_row, B);
  arma::colvec w_T_temp(B);
  arma::colvec t_T(n_row);
  arma::colvec p_temp;
  arma::colvec q;
  arma::mat Y_pred_temp;
  
  for (int i = 0; i < PLS_term_selected; ++i) {
    for (int j = 0; j < B; ++j) {
      X_matrix_temp = newdata.cols(sum(X_dim[Range(0, j)]), sum(X_dim[Range(0, j+1)]) - 1);
      w_temp = as<arma::mat>(x_weight[j]).col(i);
      t_temp = X_matrix_temp*w_temp/(sqrt(X_dim[j+1]));
      t_cbind.col(j) = t_temp;
    }
    w_T_temp = x_super_weight.col(i);
    t_T = t_cbind*w_T_temp;
    t_T_all.col(i) = t_T;
    for (int j = 0; j < B; ++j) {
      X_matrix_temp = newdata.cols(sum(X_dim[Range(0, j)]), sum(X_dim[Range(0, j+1)]) - 1);
      p_temp = as<arma::mat>(x_loading[j]).col(i);
      newdata.cols(sum(X_dim[Range(0, j)]), sum(X_dim[Range(0, j+1)]) - 1) = X_matrix_temp - t_T*p_temp.t();
    }
    q = y_weight.col(i);
    Y_pred_temp = t_T*q.t();
    Y_pred = Y_pred + Y_pred_temp;
  }
  
  arma::mat Y_pred_final(n_row, F_col);
  if (center[0]) {
    if (scale[0]) {
      // center = 1 and scale = 1
      arma::mat temp_predict = Y_pred*Y_col_sd + Y_col_mean;
      Y_pred_final = temp_predict;
    } else {
      // center = 1 and scale = 0
      arma::mat temp_predict = Y_pred + Y_col_mean;
      Y_pred_final = temp_predict;
    }
  } else {
    if (scale[0]) {
      // center = 0 and scale = 1
      arma::mat temp_predict = Y_pred*Y_col_sd;
      Y_pred_final = temp_predict;
    }
  } // center = 0 and scale = 1, no action 
  
  List output = List::create(_["Y_pred"] = Y_pred_final,
                             _["NewX_super_score"] = t_T_all);
  return (output);
}


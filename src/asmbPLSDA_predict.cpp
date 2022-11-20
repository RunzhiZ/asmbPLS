#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
List asmbPLSDA_predict(List asmbPLSDA_results,
                       arma::mat newdata, 
                       int PLS_term_selected, 
                       Nullable<String> Method = R_NilValue) {
  
  Function Euclidean_distance = Environment::namespace_env("asmbPLS")["Euclidean_distance"];
  Function Mahalanobis_distance = Environment::namespace_env("asmbPLS")["Mahalanobis_distance"];
  Function PCA_Mahalanobis_distance = Environment::namespace_env("asmbPLS")["PCA_Mahalanobis_distance"];
  
  // define
  NumericVector X_dim = asmbPLSDA_results["X_dim"];
  String outcome_type = asmbPLSDA_results["Outcome_type"];
  LogicalVector if_center = asmbPLSDA_results["center"];
  LogicalVector if_scale = asmbPLSDA_results["scale"];
  arma::rowvec X_col_mean = as<arma::rowvec>(asmbPLSDA_results["X_col_mean"]);
  arma::rowvec X_col_sd = as<arma::rowvec>(asmbPLSDA_results["X_col_sd"]);
  arma::mat x_super_score = as<arma::mat>(asmbPLSDA_results["x_super_score"]);
  arma::mat F_matrix = as<arma::mat>(asmbPLSDA_results["Y_group"]);
  x_super_score = x_super_score.cols(0, PLS_term_selected - 1);
  
  // newdata scale
  int E_col = newdata.n_cols;
  if (if_center[0]) {
    if (if_scale[0]) {
      // if_center = 1 and if_scale = 1
      for (int i = 0; i < E_col; ++i) {
        arma::mat temp = (newdata.col(i) - arma::as_scalar(X_col_mean.col(i)))/arma::as_scalar(X_col_sd.col(i));
        newdata.col(i) = temp;
      }
    } else {
      // if_center = 1 and if_scale = 0
      for (int i = 0; i < E_col; ++i) {
        arma::mat temp = (newdata.col(i) - arma::as_scalar(X_col_mean.col(i)));
        newdata.col(i) = temp;
      }
    }
  } else {
    if (if_scale[0]) {
      // if_center = 0 and if_scale = 1
      for (int i = 0; i < E_col; ++i) {
        arma::mat temp = newdata.col(i)/arma::as_scalar(X_col_sd.col(i));
        newdata.col(i) = temp;
      }
    }
  } // if_center = 0 and if_scale = 1, no action 
  
  // Method used
  String Method_used;
  if (Method.isNotNull()) {
    String Method_temp(Method);
    Method_used = Method_temp;
  } else {
    if (outcome_type == "binary") {
      Method_used = "fixed_cutoff";
    } else {
      Method_used = "Max_Y";
    }
  }
  
  // variables for convenient
  int B = X_dim.length();
  X_dim.push_front(0);
  int F_col = F_matrix.n_cols;
  int n_row = newdata.n_rows;
  int n_row_fit = x_super_score.n_rows;
  arma::mat y_weight = as<arma::mat>(asmbPLSDA_results["y_weight"]);
  arma::mat Y_pred(n_row, F_col, arma::fill::zeros);
  arma::mat Y_fit(n_row_fit, F_col, arma::fill::zeros);
  arma::mat t_T_all(n_row, PLS_term_selected);
  
  // temp variables
  arma::mat X_matrix_temp;
  arma::colvec w_temp;
  arma::colvec t_temp(n_row);
  arma::mat t_cbind(n_row, B);
  arma::colvec w_T_temp(B);
  arma::colvec t_T(n_row);
  arma::colvec p_temp;
  arma::colvec q;
  arma::mat Y_pred_temp;
  arma::mat Y_fit_temp;
  
  for (int i = 0; i < PLS_term_selected; ++i) {
    for (int j = 0; j < B; ++j) {
      X_matrix_temp = newdata.cols(sum(X_dim[Range(0, j)]), sum(X_dim[Range(0, j+1)]) - 1);
      w_temp = as<arma::mat>(as<List>(asmbPLSDA_results["x_weight"])[j]).col(i);
      t_temp = X_matrix_temp*w_temp/(sqrt(X_dim[j+1]));
      t_cbind.col(j) = t_temp;
    }
    w_T_temp = as<arma::mat>(asmbPLSDA_results["x_super_weight"]).col(i);
    t_T = t_cbind*w_T_temp;
    t_T_all.col(i) = t_T;
    for (int j = 0; j < B; ++j) {
      X_matrix_temp = newdata.cols(sum(X_dim[Range(0, j)]), sum(X_dim[Range(0, j+1)]) - 1);
      p_temp = as<arma::mat>(as<List>(asmbPLSDA_results["x_loading"])[j]).col(i);
      newdata.cols(sum(X_dim[Range(0, j)]), sum(X_dim[Range(0, j+1)]) - 1) = X_matrix_temp - t_T*p_temp.t();
    }
    q = y_weight.col(i);
    arma::colvec x_super_score_fit_temp = x_super_score.col(i);
    Y_fit_temp = x_super_score_fit_temp*q.t();
    Y_fit = Y_fit + Y_fit_temp;
    Y_pred_temp = t_T*q.t();
    Y_pred = Y_pred + Y_pred_temp;
  }
  
  arma::mat Y_pred_output(n_row, F_col);
  // Outcome: binary
  if (outcome_type == "binary") {
    // Fixed cutoff
    if (Method_used == "fixed_cutoff") {
      arma::mat temp(n_row, F_col);
      arma::colvec F_group = unique(F_matrix);
      double fixed_cutoff = mean(F_group);
      for (int i = 0; i < n_row; ++i) {
        if(Y_pred(i, 0) > fixed_cutoff) {
          temp(i, 0) = F_group[1];
        } else {
          temp(i, 0) = F_group[0];
        }
      }
      Y_pred_output = temp;
    }
    // Euclidean distance on X super score
    if (Method_used == "Euclidean_distance_X") {
      Y_pred_output = as<arma::mat>(Euclidean_distance(x_super_score, t_T_all, F_matrix, outcome_type));
    }
    // Mahalanobis distance on X super score
    if (Method_used == "Mahalanobis_distance_X") {
      Y_pred_output = as<arma::mat>(Mahalanobis_distance(x_super_score, t_T_all, F_matrix, outcome_type));
    }
  }
  
  // Outcome: more than 2 levels
  if (outcome_type == "morethan2levels") {
    // Max Y
    if (Method_used == "Max_Y") {
      for (int i = 0; i < n_row; ++i) {
        arma::rowvec row_temp = Y_pred.row(i);
        int index_max = arma::index_max(row_temp);
        Y_pred_output(i, index_max) = 1;
      }
    }
    // Euclidean distance on X super score
    if (Method_used == "Euclidean_distance_X") {
      Y_pred_output = as<arma::mat>(Euclidean_distance(x_super_score, t_T_all, F_matrix, outcome_type));
    }
    // Mahalanobis distance on X super score
    if (Method_used == "Mahalanobis_distance_X") {
      Y_pred_output = as<arma::mat>(Mahalanobis_distance(x_super_score, t_T_all, F_matrix, outcome_type));
    }
    // Euclidean distance on Y
    if (Method_used == "Euclidean_distance_Y") {
      Y_pred_output = as<arma::mat>(Euclidean_distance(Y_fit, Y_pred, F_matrix, outcome_type));
    }
    // PCA + Mahalanobis distance on Y
    if (Method_used == "PCA_Mahalanobis_distance_Y") {
      Y_pred_output = as<arma::mat>(PCA_Mahalanobis_distance(Y_fit, Y_pred));
    }
    
  }
  
  List output = List::create(_["Y_pred"] = Y_pred_output,
                             _["Y_pred_numeric"] = Y_pred,
                             _["NewX_super_score"] = t_T_all,
                             _["Method"] = Method_used);
  
  return (output);
}


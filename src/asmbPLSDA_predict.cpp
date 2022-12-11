#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
List asmbPLSDA_predict(List asmbPLSDA_results,
                       arma::mat newdata, 
                       int PLS_term_selected, 
                       String Method) {
  
  Function Euclidean_distance = Environment::namespace_env("asmbPLS")["Euclidean_distance"];
  Function Mahalanobis_distance = Environment::namespace_env("asmbPLS")["Mahalanobis_distance"];
  Function PCA_Mahalanobis_distance = Environment::namespace_env("asmbPLS")["PCA_Mahalanobis_distance"];
  
  // Extract information from asmbPLS-DA fit results
  NumericVector X_dim = asmbPLSDA_results["X_dim"];
  String outcome_type = asmbPLSDA_results["Outcome_type"];
  LogicalVector center = asmbPLSDA_results["center"];
  LogicalVector scale = asmbPLSDA_results["scale"];
  arma::rowvec X_col_mean = as<arma::rowvec>(asmbPLSDA_results["X_col_mean"]);
  arma::rowvec X_col_sd = as<arma::rowvec>(asmbPLSDA_results["X_col_sd"]);
  arma::rowvec Y_col_mean = as<arma::rowvec>(asmbPLSDA_results["Y_col_mean"]);
  List x_weight = as<List>(asmbPLSDA_results["X_weight"]);
  arma::mat x_super_weight = as<arma::mat>(asmbPLSDA_results["X_super_weight"]);
  arma::mat x_super_score = as<arma::mat>(asmbPLSDA_results["X_super_score"]);
  List x_loading = as<List>(asmbPLSDA_results["X_loading"]);
  arma::mat F_matrix = as<arma::mat>(asmbPLSDA_results["Y_group"]);
  x_super_score = x_super_score.cols(0, PLS_term_selected - 1);
  arma::mat y_weight = as<arma::mat>(asmbPLSDA_results["Y_weight"]);
  
  int E_col = newdata.n_cols; // Number of features
  int B = X_dim.length(); // Number of blocks
  int n_row_fit = x_super_score.n_rows; // Number of samples for fit data
  int F_col = F_matrix.n_cols; // Number of groups of samples
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
  arma::mat Y_fit(n_row_fit, F_col, arma::fill::zeros); // Prediction for fit data
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
  arma::mat Y_fit_temp;
  
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
    Y_fit_temp = x_super_score.col(i)*q.t();
    Y_fit = Y_fit + Y_fit_temp;
    Y_pred_temp = t_T*q.t();
    Y_pred = Y_pred + Y_pred_temp;
  }
  
  // center = 1, add the mean of Y back
  if (center[0]) {
    for (int i = 0; i < F_col; ++i) {
      arma::mat temp_fit = Y_fit.col(i) + arma::as_scalar(Y_col_mean(i));
      arma::mat temp_predict = Y_pred.col(i) + arma::as_scalar(Y_col_mean(i));
      Y_fit.col(i) = temp_fit;
      Y_pred.col(i) = temp_predict;
    }
  }
  
  arma::mat Y_pred_output(n_row, F_col, arma::fill::zeros);
  // Outcome: binary
  if (outcome_type == "binary") {
    // Fixed cutoff
    if (Method == "fixed_cutoff") {
      arma::uvec ids = find(Y_pred > 0.5);
      Y_pred_output.rows(ids).fill(1);
    }
    // Euclidean distance on X super score
    if (Method == "Euclidean_distance_X") {
      Y_pred_output = as<arma::mat>(Euclidean_distance(x_super_score, t_T_all, F_matrix, outcome_type));
    }
    // Mahalanobis distance on X super score
    if (Method == "Mahalanobis_distance_X") {
      Y_pred_output = as<arma::mat>(Mahalanobis_distance(x_super_score, t_T_all, F_matrix, outcome_type));
    }
  }
  
  // Outcome: more than 2 levels
  if (outcome_type == "multiclass") {
    // Max Y
    if (Method == "Max_Y") {
      for (int i = 0; i < n_row; ++i) {
        arma::rowvec row_temp = Y_pred.row(i);
        int index_max = arma::index_max(row_temp);
        Y_pred_output(i, index_max) = 1;
      }
    }
    // Euclidean distance on X super score
    if (Method == "Euclidean_distance_X") {
      Y_pred_output = as<arma::mat>(Euclidean_distance(x_super_score, t_T_all, F_matrix, outcome_type));
    }
    // Mahalanobis distance on X super score
    if (Method == "Mahalanobis_distance_X") {
      Y_pred_output = as<arma::mat>(Mahalanobis_distance(x_super_score, t_T_all, F_matrix, outcome_type));
    }
    // Euclidean distance on Y
    if (Method == "Euclidean_distance_Y") {
      Y_pred_output = as<arma::mat>(Euclidean_distance(Y_fit, Y_pred, F_matrix, outcome_type));
    }
    // PCA + Mahalanobis distance on Y
    if (Method == "PCA_Mahalanobis_distance_Y") {
      Y_pred_output = as<arma::mat>(PCA_Mahalanobis_distance(Y_fit, Y_pred));
    }
    
  }
  
  List output = List::create(_["Y_pred"] = Y_pred_output,
                             _["Y_pred_numeric"] = Y_pred,
                             _["NewX_super_score"] = t_T_all,
                             _["Method"] = Method);
  
  return (output);
}


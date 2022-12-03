#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
List asmbPLSDA_binary_fit(arma::mat E_matrix, 
                          arma::mat F_matrix, 
                          int PLS_term, 
                          NumericVector X_dim, 
                          arma::mat percent,
                          LogicalVector center, 
                          LogicalVector scale,
                          int maxiter) {
  
  Environment stats("package:stats");
  Function quantile_f = stats["quantile"];
  Function weight_sparse = Environment::namespace_env("asmbPLS")["weight_sparse"];
  
  int E_col = E_matrix.n_cols;
  arma::mat F_matrix_origin = F_matrix;
  // Different groups of samples
  arma::colvec F_group = unique(F_matrix);
  arma::uvec g1_index = find(F_matrix == F_group[0]);
  arma::uvec g2_index = find(F_matrix == F_group[1]);
  arma::mat E_matrix_1 = E_matrix.rows(g1_index);
  arma::mat E_matrix_2 = E_matrix.rows(g2_index);
  
  // column mean and sd
  arma::rowvec col_mean;
  arma::rowvec col_sd;
  double Y_mean;
  String outcome_type = "binary";
  
  //Center and scale
  if(center[0]) {
    // center for X
    arma::rowvec col_mean_1 = mean(E_matrix_1, 0);
    arma::rowvec col_mean_2 = mean(E_matrix_2, 0);
    col_mean = (col_mean_1 + col_mean_2)/2;
    for (int i = 0; i < E_col; ++i) {
      arma::mat temp = E_matrix.col(i) - arma::as_scalar(col_mean.col(i));
      E_matrix.col(i) = temp;
    }
    // center for Y
    Y_mean = arma::as_scalar(mean(F_group, 0));
    arma::mat temp = F_matrix - Y_mean;
    F_matrix = temp;
  }
  if(scale[0]) {
    col_sd = stddev(E_matrix, 0, 0);
    col_sd.elem(find(col_sd == 0)).ones();
    for (int i = 0; i < E_col; ++i) {
      arma::mat temp = E_matrix.col(i)/arma::as_scalar(col_sd.col(i));
      E_matrix.col(i) = temp;
    }
  }
  
  //scaled data
  arma::mat E_matrix_scaled = E_matrix;
  arma::mat F_matrix_scaled = F_matrix;
  
  // variables for convenient
  int B = X_dim.length();
  NumericVector X_dim_o = X_dim;
  X_dim.push_front(0);
  int n_row = E_matrix.n_rows;
  int F_col = F_matrix.n_cols;
  
  // variables for output
  List x_weight(B);
  List x_score(B);
  List x_loading(B);
  
  arma::mat x_weight_mat(sum(X_dim), PLS_term);
  arma::mat x_score_mat(n_row*B, PLS_term);
  arma::mat x_loading_mat(sum(X_dim), PLS_term);
  
  arma::mat x_super_weight(B, PLS_term);
  arma::mat x_super_score(n_row, PLS_term);
  arma::rowvec y_weight(PLS_term);
  arma::mat y_score(n_row, PLS_term);
  
  // temp variables
  arma::mat X_matrix_temp;
  
  arma::colvec w_temp;
  arma::colvec t_temp(n_row);
  double l_temp;
  arma::mat t_cbind(n_row, B);
  arma::colvec w_T_temp(B);
  arma::colvec p_temp;
  arma::colvec q(F_col);
  arma::colvec u(n_row);
  arma::colvec t_T(n_row);
  
  // asmbPLS fit
  arma::mat t_diff = arma::ones<arma::colvec>(n_row);
  
  for (int i = 0; i < PLS_term; ++i) {
    u = F_matrix;
    t_T = F_matrix;
    int n_iter = 0;
    while (accu(abs(t_diff)) > 0.00001 && n_iter <= maxiter) {
      
      for (int j = 0; j < B; ++j) {
        X_matrix_temp = E_matrix.cols(sum(X_dim[Range(0, j)]), sum(X_dim[Range(0, j+1)]) - 1);
        w_temp = X_matrix_temp.t()*u/arma::as_scalar(u.t()*u);
        if (percent(i, j) == 0) {
          l_temp = 0;
        } else {
          l_temp = as<double>(quantile_f(abs(w_temp), percent(i, j)));
        }
        w_temp = as<arma::mat>(weight_sparse(w_temp, l_temp));
        w_temp = w_temp/sqrt(arma::as_scalar(w_temp.t()*w_temp));
        t_temp = X_matrix_temp*w_temp/sqrt(X_dim[j+1]);
        t_cbind.col(j) = t_temp;
        
        //save x_weight, x_score
        x_weight_mat.submat(sum(X_dim[Range(0, j)]), i, sum(X_dim[Range(0, j+1)]) - 1, i) = w_temp;
        x_score_mat.submat((j)*n_row, i, (j+1)*n_row-1, i) = t_temp;
      }
      w_T_temp = t_cbind.t()*u/arma::as_scalar(u.t()*u);
      w_T_temp = w_T_temp/sqrt(arma::as_scalar(w_T_temp.t()*w_T_temp));
      t_diff = t_T - t_cbind*w_T_temp/arma::as_scalar(w_T_temp.t()*w_T_temp);
      t_T = t_cbind*w_T_temp/arma::as_scalar(w_T_temp.t()*w_T_temp);
      q = F_matrix.t()*t_T/arma::as_scalar(t_T.t()*t_T);
      u = F_matrix*q/arma::as_scalar(q.t()*q);
      // save matrix
      x_super_weight.col(i) = w_T_temp;
      x_super_score.col(i) = t_T;
      y_weight.col(i) = q;
      y_score.col(i) = u;
      // number of iterations
      n_iter = n_iter + 1;
    } 
    // make while function run in the next loop
    t_diff = arma::ones<arma::mat>(n_row, 1);
    
    // X and Y deflation
    for (int j = 0; j < B; ++j) {
      X_matrix_temp = E_matrix.cols(sum(X_dim[Range(0, j)]), sum(X_dim[Range(0, j+1)]) - 1);
      p_temp = X_matrix_temp.t()*t_T/arma::as_scalar(t_T.t()*t_T);
      X_matrix_temp = X_matrix_temp - t_T*(p_temp.t());
      E_matrix.cols(sum(X_dim[Range(0, j)]), sum(X_dim[Range(0, j+1)]) - 1) = X_matrix_temp;
      // save x_loading
      x_loading_mat.submat(sum(X_dim[Range(0, j)]), i, sum(X_dim[Range(0, j+1)]) - 1, i) = p_temp;
    }
    F_matrix = F_matrix - t_T*(q.t());
  }
  
  for (int j = 0; j < B; ++j) {
    x_weight[j] = x_weight_mat.rows(sum(X_dim[Range(0, j)]), sum(X_dim[Range(0, j+1)]) - 1);
    x_score[j] = x_score_mat.rows(j*n_row, (j+1)*n_row-1);
    x_loading[j] = x_loading_mat.rows(sum(X_dim[Range(0, j)]), sum(X_dim[Range(0, j+1)]) - 1);
  }
  List output = List::create(_["X_dim"] = X_dim_o,
                             _["x_weight"] = x_weight,
                             _["x_score"] = x_score,
                             _["x_loading"] = x_loading,
                             _["x_super_weight"] = x_super_weight,
                             _["x_super_score"] = x_super_score,
                             _["y_weight"] = y_weight,
                             _["y_score"] = y_score,
                             _["X_scaled"] = E_matrix_scaled,
                             _["Y_scaled"] = F_matrix_scaled,
                             _["X_col_mean"] = col_mean,
                             _["Y_col_mean"] = Y_mean,
                             _["X_col_sd"] = col_sd,
                             _["center"] = center,
                             _["scale"] = scale,
                             _["Outcome_type"] = outcome_type,
                             _["Y_group"] = F_matrix_origin);
  
  return (output);
}
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
List asmbPLSDA_CV(arma::mat E_matrix, 
                  arma::mat F_matrix, 
                  int PLS_term, 
                  NumericVector X_dim, 
                  arma::mat quantile_table,
                  String outcome_type,
                  String Method,
                  int K,
                  int ncv,
                  Nullable<LogicalVector> center = R_NilValue, 
                  Nullable<LogicalVector> scale = R_NilValue) {
  
  Function asmbPLSDA_fit = Environment::namespace_env("asmbPLS")["asmbPLSDA_fit"];
  Function asmbPLSDA_predict = Environment::namespace_env("asmbPLS")["asmbPLSDA_predict"];
  Function CV_index_binary = Environment::namespace_env("asmbPLS")["CV_index_binary"];
  Function CV_index_morethan2levels = Environment::namespace_env("asmbPLS")["CV_index_morethan2levels"];
  Function Results_comparison_accuracy = Environment::namespace_env("asmbPLS")["Results_comparison_accuracy"];
  
  // check if center or scale
  LogicalVector if_center;
  LogicalVector if_scale;
  if(center.isNotNull()){
    LogicalVector center_temp(center);
    if_center = center_temp[0];
  } else {
    if_center = true;
  }
  if(scale.isNotNull()){
    LogicalVector scale_temp(scale);
    if_scale = scale_temp[0];
  } else {
    if_scale = true;
  }
  
  arma::mat quantile_optimal_table(PLS_term, X_dim.size());
  int n_quantile_comb = quantile_table.n_rows;
  arma::mat quantile_table_CV(PLS_term, X_dim.size() + 1);
  arma::mat results_CV(n_quantile_comb, K);
  List CV_results(PLS_term);
  arma::uvec validation_index;
  arma::uvec training_index;
  
  List CV_index_results(ncv);
  
  arma::mat quantile_table_accuracy(n_quantile_comb, X_dim.size() + 1);
  quantile_table_accuracy.cols(0, X_dim.size() - 1) = quantile_table;
  
  if (outcome_type == "binary") {
    for (int n = 0; n < ncv; ++n) {
      CV_index_results[n] = CV_index_binary(F_matrix, K);
    }
  }
  
  if (outcome_type == "morethan2levels") {
    for (int n = 0; n < ncv; ++n) {
      CV_index_results[n] = CV_index_morethan2levels(F_matrix, K);
    }
  }
  
  for (int i = 0; i < PLS_term; ++i) {
    
    arma::mat results_CV_n(n_quantile_comb, ncv);
    
    for (int n = 0; n < ncv; ++n) {
      List CV_index_results_temp = CV_index_results[n];
      for (int j = 0; j < K; ++j) {
        // For different fold
        List CV_index_temp = CV_index_results_temp(j);
        validation_index = as<arma::uvec>(CV_index_temp["validation_index"]);
        training_index = as<arma::uvec>(CV_index_temp["training_index"]);  
        
        
        // obtain validation and training sets
        arma::mat E_matrix_validation = E_matrix.rows(validation_index);
        arma::mat F_matrix_validation = F_matrix.rows(validation_index);
        arma::mat E_matrix_training = E_matrix.rows(training_index);
        arma::mat F_matrix_training = F_matrix.rows(training_index);
        
        for (int l = 0; l < n_quantile_comb; ++l) {
          quantile_table_accuracy.submat(l, 0, l, X_dim.size() - 1) = quantile_table.row(l);
          quantile_optimal_table.row(i) = quantile_table.row(l);
          arma::mat quantile_temp =  quantile_optimal_table.rows(0, i);
          // fit model using training set
          List asmbPLSDA_fit_results = asmbPLSDA_fit(E_matrix_training, F_matrix_training, i+1, X_dim, quantile_temp, outcome_type, if_center, if_scale);
          List asmbPLSDA_predict_results = asmbPLSDA_predict(E_matrix_validation, i+1, asmbPLSDA_fit_results, Method);
          arma::mat Y_pred = asmbPLSDA_predict_results["Y_pred"];
          double accuracy = as<double>(Results_comparison_accuracy(Y_pred, F_matrix_validation));
          results_CV(l, j) = accuracy;
        }
      }
      results_CV_n.col(n) = mean(results_CV, 1);
    }
    arma::colvec results_CV_mean = mean(results_CV_n, 1);
    quantile_table_accuracy.col(X_dim.size()) = results_CV_mean;
    int index_max_accuracy = results_CV_mean.index_max();
    quantile_optimal_table.row(i) = quantile_table.row(index_max_accuracy);
    quantile_table_CV.submat(i, 0, i, X_dim.size() - 1) = quantile_table.row(index_max_accuracy);
    quantile_table_CV(i, X_dim.size()) = results_CV_mean(index_max_accuracy);
    
    CV_results(i) = quantile_table_accuracy;
  }
  
  arma::colvec PLS_accuracy = quantile_table_CV.col(X_dim.size());
  int optimal_nPLS = 1;
  if(PLS_term > 1) {
    for (int i = 0; i < PLS_term - 1; ++i) {
      double current_accuracy = arma::as_scalar(PLS_accuracy.row(i));
      double next_accuracy = arma::as_scalar(PLS_accuracy.row(i + 1));
      if(next_accuracy > current_accuracy + 0.01) {
        optimal_nPLS = optimal_nPLS + 1;
      } else {break;}
    }
  }
  
  
  List output = List::create(_["quantile_table_CV"] = quantile_table_CV,
                             _["optimal_nPLS"] = optimal_nPLS,
                             _["CV_results"] = CV_results);
  
  return(output);
}
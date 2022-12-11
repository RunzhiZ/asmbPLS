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
                  String Measure,
                  int K,
                  int ncv,
                  double expected_measure_increase,
                  LogicalVector center, 
                  LogicalVector scale,
                  int maxiter) {
  
  Function asmbPLSDA_fit = Environment::namespace_env("asmbPLS")["asmbPLSDA_fit"];
  Function asmbPLSDA_predict = Environment::namespace_env("asmbPLS")["asmbPLSDA_predict"];
  Function CV_index_binary = Environment::namespace_env("asmbPLS")["CV_index_binary"];
  Function CV_index_multiclass = Environment::namespace_env("asmbPLS")["CV_index_multiclass"];
  Function Results_comparison_measure = Environment::namespace_env("asmbPLS")["Results_comparison_measure"];
  
  double n_quantile_comb = quantile_table.n_rows; // Number of quantile combination used for CV
  double n_group = F_matrix.n_cols; // Number of groups of samples
  arma::mat quantile_table_CV(PLS_term, X_dim.size() + 5); // Table to save the best quantile combination and the corresponding measures
  
  List CV_index_results(ncv);
  // obtain the index used for CV
  if (outcome_type == "binary") {
    for (int n = 0; n < ncv; ++n) {
      CV_index_results[n] = CV_index_binary(F_matrix, K);
    }
  }
  if (outcome_type == "multiclass") {
    for (int n = 0; n < ncv; ++n) {
      CV_index_results[n] = CV_index_multiclass(F_matrix, K);
    }
  }
  
  double measure_selected;
  if(Measure == "accuracy") {
    measure_selected = 0;
  }
  if(Measure == "B_accuracy") {
    measure_selected = 1;
  }
  if(Measure == "precision") {
    measure_selected = 2;
  }
  if(Measure == "recall") {
    measure_selected = 3;
  }
  if(Measure == "F1") {
    measure_selected = 4;
  }
  
  for (int i = 0; i < PLS_term; ++i) {
    arma::mat results_CV_summary_n(n_quantile_comb, ncv, arma::fill::zeros);
    for (int n = 0; n < ncv; ++n) {
      List CV_index_results_temp = CV_index_results[n];
      // K folds
      for (int j = 0; j < K; ++j) {
        List CV_index_temp = CV_index_results_temp(j);
        arma::uvec validation_index = as<arma::uvec>(CV_index_temp["validation_index"]);
        arma::uvec training_index = as<arma::uvec>(CV_index_temp["training_index"]);  
        
        // obtain validation and training sets
        arma::mat E_matrix_validation = E_matrix.rows(validation_index);
        arma::mat F_matrix_validation = F_matrix.rows(validation_index);
        arma::mat E_matrix_training = E_matrix.rows(training_index);
        arma::mat F_matrix_training = F_matrix.rows(training_index);
        
        //calculate overall/balanced accuracy using different quantile combinations
        for (int l = 0; l < n_quantile_comb; ++l) {
          quantile_table_CV.submat(i, 0, i, X_dim.size() - 1) = quantile_table.row(l);
          arma::mat quantile_temp =  quantile_table_CV.submat(0, 0, i, X_dim.size() - 1);
          // fit model using training set
          List asmbPLSDA_fit_results = asmbPLSDA_fit(E_matrix_training, F_matrix_training, i+1, X_dim, quantile_temp, outcome_type, center, scale, maxiter);
          List asmbPLSDA_predict_results = asmbPLSDA_predict(asmbPLSDA_fit_results, E_matrix_validation, i+1, Method);
          arma::mat Y_pred = asmbPLSDA_predict_results["Y_pred"];
          arma::rowvec measure = as<arma::rowvec>(Results_comparison_measure(Y_pred, F_matrix_validation, outcome_type));
          
          double temp = results_CV_summary_n(l, n);
          temp = temp + measure(measure_selected);
          results_CV_summary_n(l, n) = temp;
        }
      }
    }
    results_CV_summary_n = results_CV_summary_n/ncv;
    //calculate the mean accuracy for n CVs
    arma::colvec results_CV_summary_mean = mean(results_CV_summary_n, 1);
    
    // find the quantile combination with the highest accuracy
    int index_max_measure = results_CV_summary_mean.index_max();
    quantile_table_CV.submat(i, 0, i, X_dim.size() - 1) = quantile_table.row(index_max_measure);
    
    arma::mat Y_pred_bind(1, n_group, arma::fill::zeros);
    arma::mat F_matrix_validation_bind(1, n_group, arma::fill::zeros);
    
    // obtain corresponding measures for the selected quantile combination
    for (int n = 0; n < ncv; ++n) {
      List CV_index_results_temp = CV_index_results[n];
      for (int j = 0; j < K; ++j) {
        // For different fold
        List CV_index_temp = CV_index_results_temp(j);
        arma::uvec validation_index = as<arma::uvec>(CV_index_temp["validation_index"]);
        arma::uvec training_index = as<arma::uvec>(CV_index_temp["training_index"]);  
        
        // obtain validation and training sets
        arma::mat E_matrix_validation = E_matrix.rows(validation_index);
        arma::mat F_matrix_validation = F_matrix.rows(validation_index);
        arma::mat E_matrix_training = E_matrix.rows(training_index);
        arma::mat F_matrix_training = F_matrix.rows(training_index);
        
        arma::mat quantile_temp =  quantile_table_CV.submat(0, 0, i, X_dim.size() - 1);
        // fit model using training set
        List asmbPLSDA_fit_results = asmbPLSDA_fit(E_matrix_training, F_matrix_training, i+1, X_dim, quantile_temp, outcome_type, center, scale, maxiter);
        List asmbPLSDA_predict_results = asmbPLSDA_predict(asmbPLSDA_fit_results, E_matrix_validation, i+1, Method);
        arma::mat Y_pred = asmbPLSDA_predict_results["Y_pred"];
        
        arma::mat temp_pred = join_vert(Y_pred_bind, Y_pred);
        arma::mat temp_validation = join_vert(F_matrix_validation_bind, F_matrix_validation);
        Y_pred_bind = temp_pred;
        F_matrix_validation_bind = temp_validation;
      }
    }
    arma::rowvec measure = as<arma::rowvec>(Results_comparison_measure(Y_pred_bind.rows(1, Y_pred_bind.n_rows - 1), F_matrix_validation_bind.rows(1, F_matrix_validation_bind.n_rows - 1), outcome_type));
    quantile_table_CV.submat(i, X_dim.size(), i, X_dim.size() + 4) = measure;
  }
  
  arma::colvec PLS_measure = quantile_table_CV.col(X_dim.size() + measure_selected);
  int optimal_nPLS = 1;
  if(PLS_term > 1) {
    for (int i = 0; i < PLS_term - 1; ++i) {
      double current_measure = arma::as_scalar(PLS_measure.row(i));
      double next_measure = arma::as_scalar(PLS_measure.row(i + 1));
      if(next_measure > current_measure + expected_measure_increase) {
        optimal_nPLS = optimal_nPLS + 1;
      } else {break;}
    }
  }
  
  List output = List::create(_["quantile_table_CV"] = quantile_table_CV,
                             _["optimal_nPLS"] = optimal_nPLS);
  
  return(output);
}
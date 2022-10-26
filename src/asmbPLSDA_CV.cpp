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
                  Nullable<int> K = R_NilValue,
                  Nullable<LogicalVector> center = R_NilValue, 
                  Nullable<LogicalVector> scale = R_NilValue,
                  Nullable<int> seed = R_NilValue) {
  
  Function asmbPLSDA_fit = Environment::namespace_env("asmbPLS")["asmbPLSDA_fit"];
  Function asmbPLSDA_predict = Environment::namespace_env("asmbPLS")["asmbPLSDA_predict"];
  Function CV_index_binary = Environment::namespace_env("asmbPLS")["CV_index_binary"];
  Function CV_index_morethan2levels = Environment::namespace_env("asmbPLS")["CV_index_morethan2levels"];
  Function Results_comparison_accuracy = Environment::namespace_env("asmbPLS")["Results_comparison_accuracy"];
  // Function asmbPLSDA_fit("asmbPLSDA_fit");
  // Function asmbPLSDA_predict("asmbPLSDA_predict");
  // Function CV_index_binary("CV_index_binary");
  // Function CV_index_morethan2levels("CV_index_morethan2levels");
  // Function Results_comparison_accuracy("Results_comparison_accuracy");
  
  
  // Set default K = 5 if no input K
  int K_input;
  if(K.isNotNull()) {
    IntegerVector K_temp(K);
    K_input = K_temp[0];
  } else{
    K_input = 5;
  }
  
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
  
  // Set default seed = 1
  int seed_input;
  if(seed.isNotNull()) {
    IntegerVector seed_temp(seed);
    seed_input = seed_temp[0];
  } else{
    seed_input = 1;
  }
  
  arma::mat quantile_optimal_table(PLS_term, X_dim.size());
  int n_quantile_comb = quantile_table.n_rows;
  arma::mat quantile_table_CV(PLS_term, X_dim.size() + 1);
  arma::mat results_CV(n_quantile_comb, K_input);
  List CV_results(PLS_term);
  arma::uvec validation_index;
  arma::uvec training_index;
  
  List CV_index_results;
  List CV_index_results_output;
  
  arma::mat quantile_table_accuracy(n_quantile_comb, X_dim.size() + 1);
  quantile_table_accuracy.cols(0, X_dim.size() - 1) = quantile_table;
  
  
  if (outcome_type == "binary") {
    CV_index_results = CV_index_binary(F_matrix, K_input, seed_input);
    CV_index_results_output = CV_index_binary(F_matrix, K_input, seed_input);
  }
  
  if (outcome_type == "morethan2levels") {
    CV_index_results = CV_index_morethan2levels(F_matrix, K_input, seed_input);
    CV_index_results_output = CV_index_morethan2levels(F_matrix, K_input, seed_input);
  }
  
  for (int i = 0; i < CV_index_results_output.size(); ++i) {
    List temp = CV_index_results_output[i];
    for (int j = 0; j < temp.size(); ++j) {
      NumericVector index_temp = temp[j];
      for (int l = 0; l < index_temp.size(); ++l) {
        index_temp[l] = index_temp[l] + 1;
      }
    }
  }
  
  for (int i = 0; i < PLS_term; ++i) {
    for (int j = 0; j < K_input; ++j) {
      // For different fold
      
      List CV_index_temp = CV_index_results(j);
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
    arma::colvec results_CV_mean = mean(results_CV, 1);
    quantile_table_accuracy.col(X_dim.size()) = results_CV_mean;
    int index_max_accuracy = results_CV_mean.index_max();
    quantile_optimal_table.row(i) = quantile_table.row(index_max_accuracy);
    quantile_table_CV.submat(i, 0, i, X_dim.size() - 1) = quantile_table.row(index_max_accuracy);
    quantile_table_CV(i, X_dim.size()) = results_CV_mean(index_max_accuracy);
    
    CV_results(i) = quantile_table_accuracy;
  }
  
  List output = List::create(_["quantile_table_CV"] = quantile_table_CV,
                             _["CV_results"] = CV_results,
                             _["CV_index"] = CV_index_results_output);
  
  return(output);
}
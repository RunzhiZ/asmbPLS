#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
List asmbPLS_CV(arma::mat E_matrix, 
                arma::mat F_matrix, 
                int PLS_term, 
                NumericVector X_dim, 
                arma::mat quantile_table,
                arma::colvec Y_indicator,
                Nullable<int> K = R_NilValue,
                Nullable<LogicalVector> only_observe = R_NilValue,
                Nullable<int> seed = R_NilValue) {
  
  Function asmbPLS_fit = Environment::namespace_env("asmbPLS")["asmbPLS_fit"];
  Function asmbPLS_predict = Environment::namespace_env("asmbPLS")["asmbPLS_predict"];
  Function CV_index = Environment::namespace_env("asmbPLS")["CV_index"];
  Function Results_comparison_MSE = Environment::namespace_env("asmbPLS")["Results_comparison_MSE"];
  // Function asmbPLS_fit("asmbPLS_fit");
  // Function asmbPLS_predict("asmbPLS_predict");
  // Function CV_index("CV_index");
  // Function Results_comparison_MSE("Results_comparison_MSE");
  
  // Set default K = 5 if no input K
  int K_input;
  if(K.isNotNull()) {
    IntegerVector K_temp(K);
    K_input = K_temp[0];
  } else{
    K_input = 5;
  }
  
  // Set default only_observe = T
  LogicalVector only_observe_input;
  if(only_observe.isNotNull()) {
    LogicalVector only_observe_temp(only_observe);
    only_observe_input = only_observe_temp[0];
  } else{
    only_observe_input = true;
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
  
  
  arma::mat quantile_table_MSE(n_quantile_comb, X_dim.size() + 1);
  quantile_table_MSE.cols(0, 1) = quantile_table;
  
  
  List CV_index_results = CV_index(Y_indicator, K_input, seed_input, only_observe_input);
  List CV_index_results_output = CV_index(Y_indicator, K_input, seed_input, only_observe_input);
  
  for (int i = 0; i < CV_index_results_output.size(); ++i) {
    List temp = CV_index_results_output[i];
    for (int j = 0; j < temp.size(); ++j) {
      NumericVector index_temp = temp[j];
      for (int l = 0; l < index_temp.size(); ++l) {
        index_temp[l] = index_temp[l] + 1;
      }
      temp[j] = index_temp;
    }
    CV_index_results_output[i] = temp;
  }
  
  for (int i = 0; i < PLS_term; ++i) {
    for (int j = 0; j < K_input; ++j) {
      // For different fold
      List CV_index_temp = CV_index_results[j];
      // obtain validation and training sets
      arma::uvec validation_index = as<arma::uvec>(CV_index_temp["validation_index"]);
      arma::uvec training_index = as<arma::uvec>(CV_index_temp["training_index"]);

      arma::mat E_matrix_validation = E_matrix.rows(validation_index);
      arma::mat F_matrix_validation = F_matrix.rows(validation_index);
      arma::mat E_matrix_training = E_matrix.rows(training_index);
      arma::mat F_matrix_training = F_matrix.rows(training_index);
      
      
      for (int l = 0; l < n_quantile_comb; ++l) {
        quantile_optimal_table.row(i) = quantile_table.row(l);
        arma::mat quantile_temp =  quantile_optimal_table.rows(0, i);
        // fit model using training set
        List asmbPLS_fit_results = asmbPLS_fit(E_matrix_training, F_matrix_training, i + 1, X_dim, quantile_temp);
        arma::mat Y_pred = as<arma::mat>(asmbPLS_predict(E_matrix_validation, i + 1, asmbPLS_fit_results));
        double MSE = as<double>(Results_comparison_MSE(Y_pred, F_matrix_validation));
        results_CV(l, j) = MSE;
      }
    }
    arma::colvec results_CV_mean = mean(results_CV, 1);
    quantile_table_MSE.col(2) = results_CV_mean;
    int index_min_MSE = results_CV_mean.index_min();
    quantile_optimal_table.row(i) = quantile_table.row(index_min_MSE);
    quantile_table_CV.submat(i, 0, i, 1) = quantile_table.row(index_min_MSE);
    quantile_table_CV(i, 2) = results_CV_mean(index_min_MSE);
    
    CV_results(i) = quantile_table_MSE;
  }
  
  List output = List::create(_["quantile_table_CV"] = quantile_table_CV,
                             _["CV_results"] = CV_results,
                             _["CV_index"] = CV_index_results_output);
  
  return(output);
}
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
List CV_index(arma::colvec Y_indicator, 
              int K_input,
              LogicalVector only_observe) {
  
  // call functions
  Function stl_sort = Environment::namespace_env("asmbPLS")["stl_sort"];
  Function sample_group = Environment::namespace_env("asmbPLS")["sample_group"];
  Environment base("package:base");
  Function sample_function = base["sample"];
  
  // Take same proportion samples from both training and validation sets
  int n_sample = Y_indicator.n_rows; // Number of samples
  arma::colvec F_group = unique(Y_indicator);
  int n_group = F_group.n_rows; // Number of groups of Y_indicator
  
  List CV_index_output(K_input);
  NumericVector index_resample_all(n_sample);
  NumericVector index_K_all(n_sample);
  NumericVector index_observe_all(n_sample);
  int n_g_start = 0;
  
  for (int i = 0; i < n_group; ++i) {
    arma::uvec index_g = find(Y_indicator == F_group[i]); // group indices
    int n_g = index_g.n_rows; // Number of sample in the specific group
    NumericVector index_resample_g = sample_function(index_g, n_g, false); // Shuffle the indices
    NumericVector index_K_g = sample_group(n_g, K_input); // Group the samples
    NumericVector index_observe(n_g);
    for (int j = 0; j < n_g; ++j) {
      index_observe[j] = F_group[i];
    }
    
    for (int j = n_g_start; j < n_g_start + index_resample_g.size(); ++j) {
      index_resample_all[j] = index_resample_g[j - n_g_start];
      index_K_all[j] = index_K_g[j - n_g_start];
      index_observe_all[j] = index_observe[j - n_g_start];
    }
    
    n_g_start = n_g_start + n_g;
  }
  
  for (int i = 0; i < K_input; ++i) {
    NumericVector validation_index;
    NumericVector training_index;
    for (int j = 0; j < index_K_all.size(); ++j) {
      if (only_observe[0]) {
        if (index_K_all[j] == i && index_observe_all[j] == 1) {
          validation_index.push_back(index_resample_all[j]);
        } 
        if (!(index_K_all[j] == i)) {
          training_index.push_back(index_resample_all[j]);
        }
      } else {
        if (index_K_all[j] == i) {
          validation_index.push_back(index_resample_all[j]);
        } else {
          training_index.push_back(index_resample_all[j]);
        }
      }
    }
    validation_index = stl_sort(validation_index);
    training_index = stl_sort(training_index);
    List output = List::create(_["validation_index"] = validation_index,
                               _["training_index"] = training_index);
    CV_index_output[i] = output;
  }
  return(CV_index_output);
}
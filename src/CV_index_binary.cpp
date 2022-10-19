#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
List CV_index_binary(arma::mat F_matrix, 
                     int K_input, 
                     int seed) {
  
  // call functions
  Function stl_sort = Environment::namespace_env("asmbPLS")["stl_sort"];
  Function sample_group = Environment::namespace_env("asmbPLS")["sample_group"];
  Environment base("package:base");
  Function sample_function = base["sample"];
  Function set_seed = base["set.seed"];
  
  // Take same proportion samples from both training and validation sets
  arma::colvec F_group = unique(F_matrix);
  int n_sample = F_matrix.n_rows;
  
  arma::uvec g1_index = find(F_matrix == F_group[0]);
  arma::uvec g2_index = find(F_matrix == F_group[1]);
  
  int n_g1 = g1_index.n_rows;
  int n_g2 = g2_index.n_rows;
  
  // set seed
  set_seed(seed);
  NumericVector g1_index_resample = sample_function(g1_index, n_g1, false);
  NumericVector g2_index_resample = sample_function(g2_index, n_g2, false);
  
  NumericVector g1_index_K = sample_group(n_g1, K_input);
  NumericVector g2_index_K = sample_group(n_g2, K_input);
  List CV_index_output(K_input);
  
  for (int i = 0; i < K_input; ++i) {
    NumericVector g1_validation_index;
    NumericVector g2_validation_index;
    for (int j = 0; j < n_g1; ++j) {
      if(g1_index_K(j) == i) {
        g1_validation_index.push_back(g1_index_resample(j));
      }
    }
    for (int l = 0; l < n_g2; ++l) {
      if(g2_index_K(l) == i) {
        g2_validation_index.push_back(g2_index_resample(l));
      }
    }
    NumericVector validation_index(g1_validation_index.size() + g2_validation_index.size());
    validation_index[seq(0, g1_validation_index.size() - 1)] = g1_validation_index;
    validation_index[seq(g1_validation_index.size(), g1_validation_index.size() + g2_validation_index.size() - 1)] = g2_validation_index;
    validation_index = stl_sort(validation_index);
    
    NumericVector all_index(n_sample);
    for (int n = 0; n < n_sample; ++n) {
      all_index[n] = n;
    }
    
    LogicalVector if_validation = in(all_index, validation_index);
    LogicalVector if_training = !(if_validation);
    
    NumericVector training_index = all_index[if_training == 1];
    List output = List::create(_["validation_index"] = validation_index,
                               _["training_index"] = training_index);
    CV_index_output[i] = output;
  }
  return(CV_index_output);
}
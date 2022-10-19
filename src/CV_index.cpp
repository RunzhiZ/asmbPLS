#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
List CV_index(arma::colvec Y_indicator, 
              int K_input, 
              int seed,
              LogicalVector only_observe) {
  
  // call functions
  Function stl_sort = Environment::namespace_env("asmbPLS")["stl_sort"];
  Function sample_group = Environment::namespace_env("asmbPLS")["sample_group"];
  Environment base("package:base");
  Function sample_function = base["sample"];
  Function set_seed = base["set.seed"];
  
  // Take same proportion samples from both training and validation sets
  int n_sample = Y_indicator.n_rows;
  List CV_index_output(K_input);
  arma::colvec F_group = unique(Y_indicator);
  
  if(F_group.n_rows == 1) {
    arma::uvec g_index = find(Y_indicator == F_group[0]);
    int n_g = g_index.n_rows;
    // set seed
    set_seed(seed);
    // resample index
    NumericVector g_index_resample = sample_function(g_index, n_g, false);
    // cut index into different groups
    NumericVector g_index_K = sample_group(n_g, K_input);
    
    for (int i = 0; i < K_input; ++i) {
      NumericVector validation_index;
      for (int j = 0; j < n_g; ++j) {
        if(g_index_K(j) == i) {
          validation_index.push_back(g_index_resample(j));
        }
      }
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
  } else {
    arma::uvec g0_index = find(Y_indicator == F_group[0]);
    arma::uvec g1_index = find(Y_indicator == F_group[1]);
    
    int n_g0 = g0_index.n_rows;
    int n_g1 = g1_index.n_rows;
    
    // set seed
    set_seed(seed);
    // resample index
    NumericVector g0_index_resample = sample_function(g0_index, n_g0, false);
    NumericVector g1_index_resample = sample_function(g1_index, n_g1, false);
    // cut index into different groups
    NumericVector g0_index_K = sample_group(n_g0, K_input);
    NumericVector g1_index_K = sample_group(n_g1, K_input);
    
    for (int i = 0; i < K_input; ++i) {
      NumericVector g0_validation_index;
      NumericVector g1_validation_index;
      for (int j = 0; j < n_g0; ++j) {
        if(g0_index_K(j) == i) {
          g0_validation_index.push_back(g0_index_resample(j));
        }
      }
      for (int l = 0; l < n_g1; ++l) {
        if(g1_index_K(l) == i) {
          g1_validation_index.push_back(g1_index_resample(l));
        }
      }
      g1_validation_index = stl_sort(g1_validation_index);
      
      NumericVector validation_index(g0_validation_index.size() + g1_validation_index.size());
      validation_index[seq(0, g0_validation_index.size() - 1)] = g0_validation_index;
      validation_index[seq(g0_validation_index.size(), g0_validation_index.size() + g1_validation_index.size() - 1)] = g1_validation_index;
      validation_index = stl_sort(validation_index);
      
      NumericVector all_index(n_sample);
      for (int n = 0; n < n_sample; ++n) {
        all_index[n] = n;
      }
      
      LogicalVector if_validation = in(all_index, validation_index);
      LogicalVector if_training = !(if_validation);
      
      NumericVector training_index = all_index[if_training == 1];
      
      if(only_observe[0]) {
        List output = List::create(_["validation_index"] = g1_validation_index,
                                   _["training_index"] = training_index);
        CV_index_output[i] = output;
      } else {
        List output = List::create(_["validation_index"] = validation_index,
                                   _["training_index"] = training_index);
        CV_index_output[i] = output;
      }
    }
  }
  return(CV_index_output);
}
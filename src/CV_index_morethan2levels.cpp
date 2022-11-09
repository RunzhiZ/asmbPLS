#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
List CV_index_morethan2levels(arma::mat F_matrix, 
                              int K_input) {
  
  // call functions
  Function stl_sort = Environment::namespace_env("asmbPLS")["stl_sort"];
  Function sample_group = Environment::namespace_env("asmbPLS")["sample_group"];
  Environment base("package:base");
  Function sample_function = base["sample"];
  
  // Take same proportion samples from both training and validation sets
  arma::colvec F_group = unique(F_matrix);
  int n_sample = F_matrix.n_rows;
  int F_col = F_matrix.n_cols;
  
  // More than 2 levels
  List g_index_resample_list(F_col);
  List g_index_K_list(F_col);
  List n_g_list(F_col);
  List CV_index_output(K_input);
  
  
  for (int i = 0; i < F_col; ++i) {
    arma::uvec g_index = find(F_matrix.col(i) == 1);
    int n_g = g_index.n_rows;
    
    NumericVector g_index_resample = sample_function(g_index, n_g, false);
    NumericVector g_index_K = sample_group(n_g, K_input);
    g_index_resample_list[i] = g_index_resample;
    g_index_K_list[i] = g_index_K;
    n_g_list[i] = n_g;
  }
  
  List g_validation_index_list(F_col);
  
  for (int i = 0; i < K_input; ++i) {
    int n_validation = 0;
    for (int j = 0; j < F_col; ++j) {
      NumericVector g_validation_index;
      int n_g = as<int>(n_g_list[j]);
      NumericVector g_index_K = g_index_K_list[j];
      NumericVector g_index_resample = g_index_resample_list[j];
      for (int l = 0; l < n_g; ++l) {
        if(g_index_K(l) == i) {
          g_validation_index.push_back(g_index_resample(l));
        }
      }
      int n_g_validation = g_validation_index.size();
      n_validation = n_validation + n_g_validation;
      g_validation_index_list(j) = g_validation_index;
    }
    
    NumericVector validation_index(n_validation);
    validation_index[seq(0, as<NumericVector>(g_validation_index_list[0]).size() - 1)] = as<NumericVector>(g_validation_index_list[0]);
    int index_n_1 = 0;
    int index_n_2 = as<NumericVector>(g_validation_index_list[0]).size();
    for (int j = 1; j < F_col; ++j) {
      index_n_1 = index_n_1 + as<NumericVector>(g_validation_index_list[j-1]).size();
      index_n_2 = index_n_2 + as<NumericVector>(g_validation_index_list[j]).size();
      validation_index[seq(index_n_1, index_n_2 - 1)] = as<NumericVector>(g_validation_index_list[j]);
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
  return(CV_index_output);
}
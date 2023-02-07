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
                int K,
                int ncv,
                LogicalVector only_observe,
                double expected_MSE_decrease,
                LogicalVector center, 
                LogicalVector scale,
                int maxiter) {
  
  Function asmbPLS_fit = Environment::namespace_env("asmbPLS")["asmbPLS_fit"];
  Function asmbPLS_predict = Environment::namespace_env("asmbPLS")["asmbPLS_predict"];
  Function CV_index = Environment::namespace_env("asmbPLS")["CV_index"];
  Function Results_comparison_MSE = Environment::namespace_env("asmbPLS")["Results_comparison_MSE"];
  
  int n_quantile_comb = quantile_table.n_rows; // Number of quantile combination used for CV
  arma::mat quantile_table_CV(PLS_term, X_dim.size() + 1); // Table to save the best quantile combination and the corresponding measures
  
  List CV_index_results(ncv);
  // obtain the index used for CV
  for (int n = 0; n < ncv; ++n) {
    CV_index_results[n] = CV_index(Y_indicator, K, only_observe);
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

        //calculate MSE using different quantile combinations
        for (int l = 0; l < n_quantile_comb; ++l) {
          quantile_table_CV.submat(i, 0, i, X_dim.size() - 1) = quantile_table.row(l);
          arma::mat quantile_temp =  quantile_table_CV.submat(0, 0, i, X_dim.size() - 1);
          // fit model using training set
          List asmbPLS_fit_results = asmbPLS_fit(E_matrix_training, F_matrix_training, i+1, X_dim, quantile_temp, center, scale, maxiter);
          List asmbPLS_predict_results = asmbPLS_predict(asmbPLS_fit_results, E_matrix_validation, i+1);
          arma::mat Y_pred = asmbPLS_predict_results["Y_pred"];
          double measure = as<double>(Results_comparison_MSE(Y_pred, F_matrix_validation));

          double temp = results_CV_summary_n(l, n);
          temp = temp + measure;
          results_CV_summary_n(l, n) = temp;
        }
      }
    }
    results_CV_summary_n = results_CV_summary_n/ncv;
    //calculate the mean accuracy for n CVs
    arma::colvec results_CV_summary_mean = mean(results_CV_summary_n, 1);

    // find the quantile combination with the lowest MSE
    int index_min_measure = results_CV_summary_mean.index_min();
    quantile_table_CV.submat(i, 0, i, X_dim.size() - 1) = quantile_table.row(index_min_measure);
    quantile_table_CV(i, X_dim.size()) = results_CV_summary_mean(index_min_measure);
  }

  arma::colvec PLS_measure = quantile_table_CV.col(X_dim.size());
  int optimal_nPLS = 1;
  if(PLS_term > 1) {
    for (int i = 0; i < PLS_term - 1; ++i) {
      double current_measure = arma::as_scalar(PLS_measure.row(i));
      double next_measure = arma::as_scalar(PLS_measure.row(i + 1));
      if(next_measure <= current_measure*(1 - expected_MSE_decrease)) {
        optimal_nPLS = optimal_nPLS + 1;
      } else {break;}
    }
  }

  List output = List::create(_["quantile_table_CV"] = quantile_table_CV,
                             _["optimal_nPLS"] = optimal_nPLS);
  return(output);
}
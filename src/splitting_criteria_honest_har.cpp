// Includes/namespaces
#include <RcppArmadillo.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
Rcpp::NumericVector splitting_criterion_honest_rcpp(const arma::vec& split_var_values,
                                                    const arma::mat& X,
                                                    const arma::vec& response,
                                                    const Rcpp::NumericVector& splits) {
  int n_splits = splits.size();
  // int n = split_var_values.size();
  NumericVector loss(n_splits);

  // Preallocate variables outside the loop to avoid repeated allocations
  arma::uvec left_indices, right_indices;
  arma::mat X_left, X_right;
  arma::vec response_left, response_right, coef_left, coef_right;
  arma::vec residuals_left, residuals_right;
  double rss_left, rss_right;

  for (int ii = 0; ii < n_splits; ++ii) {
    double split_value = splits[ii];

    // Use Armadillo\'s efficient subsetting mechanism (find function)
    arma::uvec left_indices = arma::find(split_var_values < split_value);
    arma::uvec right_indices = arma::find(split_var_values >= split_value);

    // Skip if either partition is empty
    if (left_indices.n_elem == 0 || right_indices.n_elem == 0) {
      loss[ii] = NA_REAL;
      continue;
    }

    // Subset the design matrix and response vector
    X_left = X.rows(left_indices);
    response_left = response.elem(left_indices);
    X_right = X.rows(right_indices);
    response_right = response.elem(right_indices);

    // Solve the linear equations directly for both partitions
    coef_left = arma::solve(X_left, response_left);
    coef_right = arma::solve(X_right, response_right);

    // Compute residuals
    residuals_left = response_left - X_left * coef_left;
    residuals_right = response_right - X_right * coef_right;

    // Calculate the sum of squared residuals using dot product
    rss_left = arma::dot(residuals_left, residuals_left);
    rss_right = arma::dot(residuals_right, residuals_right);

    // Store the total loss for this split
    loss[ii] = rss_left + rss_right;
  }

  return loss;
}
//
// NumericVector splitting_criterion_honest_rcpp(NumericVector split_var_values,
//                                               arma::mat X,
//                                               arma::vec response,
//                                               NumericVector splits) {
//   int n_splits = splits.size();
//   int n = split_var_values.size();
//   NumericVector loss(n_splits);
//
//   for (int ii = 0; ii < n_splits; ii++) {
//     double split_value = splits[ii];
//
//     // Vectors to hold left and right indices
//     std::vector<int> left_indices, right_indices;
//
//     // Partition the data based on the split value
//     for (int i = 0; i < n; i++) {
//       if (split_var_values[i] < split_value) {
//         left_indices.push_back(i);
//       } else {
//         right_indices.push_back(i);
//       }
//     }
//
//     // If either partition is empty, skip this split
//     if (left_indices.size() == 0 || right_indices.size() == 0) {
//       loss[ii] = NA_REAL;
//       continue;
//     }
//
//     // Extract left and right subsets of X and response
//     arma::mat X_left(left_indices.size(), X.n_cols);
//     arma::vec response_left(left_indices.size());
//     arma::mat X_right(right_indices.size(), X.n_cols);
//     arma::vec response_right(right_indices.size());
//
//     for (size_t i = 0; i < left_indices.size(); i++) {
//       X_left.row(i) = X.row(left_indices[i]);
//       response_left[i] = response[left_indices[i]];
//     }
//     for (size_t i = 0; i < right_indices.size(); i++) {
//       X_right.row(i) = X.row(right_indices[i]);
//       response_right[i] = response[right_indices[i]];
//     }
//
//     // Fit the models using QR decomposition (equivalent to .lm.fit in R)
//     arma::vec coef_left = arma::solve(X_left, response_left);
//     arma::vec coef_right = arma::solve(X_right, response_right);
//
//     // Calculate residuals
//     arma::vec residuals_left = response_left - X_left * coef_left;
//     arma::vec residuals_right = response_right - X_right * coef_right;
//
//     // Compute the sum of squared residuals
//     double rss_left = arma::dot(residuals_left, residuals_left);
//     double rss_right = arma::dot(residuals_right, residuals_right);
//
//     // Store the sum of residuals squared for the split
//     loss[ii] = rss_left + rss_right;
//   }
//
//   return loss;
// }

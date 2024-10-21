//Includes/namespaces
#include <RcppArmadillo.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
List splitting_criterion_honest_intercept_only_cpp(NumericVector split_var_values, NumericVector response, NumericVector splits) {
  int n_splits = splits.size();
  int n = split_var_values.size();
  NumericVector loss(n_splits);

  // Iterate over each split point
  for (int ii = 0; ii < n_splits; ii++) {
    double split_value = splits[ii];

    // Variables for calculating means and losses
    double sum_left = 0.0, sum_right = 0.0;
    int count_left = 0, count_right = 0;

    // Calculate sums for left and right partitions
    for (int i = 0; i < n; i++) {
      if (split_var_values[i] < split_value) {
        sum_left += response[i];
        count_left++;
      } else {
        sum_right += response[i];
        count_right++;
      }
    }

    // Calculate means for left and right partitions
    double mean_left = count_left > 0 ? sum_left / count_left : 0.0;
    double mean_right = count_right > 0 ? sum_right / count_right : 0.0;

    // Calculate SSE for left and right partitions
    double loss_left = 0.0, loss_right = 0.0;
    for (int i = 0; i < n; i++) {
      if (split_var_values[i] < split_value) {
        loss_left += pow(response[i] - mean_left, 2);
      } else {
        loss_right += pow(response[i] - mean_right, 2);
      }
    }

    // Store the total loss for this split point
    loss[ii] = loss_left + loss_right;
  }

  // Find the index of the minimum loss
  double min_loss = loss[0];
  int min_index = 0;
  for (int ii = 1; ii < n_splits; ii++) {
    if (loss[ii] < min_loss) {
      min_loss = loss[ii];
      min_index = ii;
    }
  }

  // Return the minimum loss and the corresponding split value
  return List::create(
    Named("sse") = min_loss,
    Named("split") = splits[min_index]
  );
}

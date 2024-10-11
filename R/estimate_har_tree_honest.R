
# Define the Rcpp function
# Load the Rcpp and RcppArmadillo libraries
library(Rcpp)
library(RcppArmadillo)

cppFunction(depends = "RcppArmadillo", code = '
Rcpp::NumericVector splitting_criterion_honest_rcpp(Rcpp::NumericVector split_var_values,
                                                    arma::mat X,
                                                    arma::vec response,
                                                    Rcpp::NumericVector splits) {
    int n_splits = splits.size();
    int n = split_var_values.size();
    Rcpp::NumericVector loss(n_splits);

    for (int ii = 0; ii < n_splits; ++ii) {
        double split_value = splits[ii];

        // Use Armadillo\'s efficient subsetting mechanism (find function)
        arma::uvec left_indices = arma::find(arma::vec(split_var_values.begin(), split_var_values.size(), false) < split_value);
        arma::uvec right_indices = arma::find(arma::vec(split_var_values.begin(), split_var_values.size(), false) >= split_value);

        // Skip if either partition is empty
        if (left_indices.n_elem == 0 || right_indices.n_elem == 0) {
            loss[ii] = NA_REAL;
            continue;
        }

        // Subset the design matrix and response vector efficiently
        arma::mat X_left = X.rows(left_indices);
        arma::vec response_left = response.elem(left_indices);
        arma::mat X_right = X.rows(right_indices);
        arma::vec response_right = response.elem(right_indices);

        // Solve the linear equations directly for both partitions
        arma::vec coef_left = arma::solve(X_left, response_left);
        arma::vec coef_right = arma::solve(X_right, response_right);

        // Compute residuals using Armadillo\'s efficient vectorized operations
        arma::vec residuals_left = response_left - X_left * coef_left;
        arma::vec residuals_right = response_right - X_right * coef_right;

        // Calculate the sum of squared residuals using dot product
        double rss_left = arma::dot(residuals_left, residuals_left);
        double rss_right = arma::dot(residuals_right, residuals_right);

        // Store the total loss for this split
        loss[ii] = rss_left + rss_right;
    }

    return loss;
}')


# # Define the Rcpp function using RcppArmadillo
# cppFunction(depends = "RcppArmadillo", '
# NumericVector splitting_criterion_honest_rcpp(NumericVector split_var_values,
#                                               arma::mat X,
#                                               arma::vec response,
#                                               NumericVector splits) {
#     int n_splits = splits.size();
#     int n = split_var_values.size();
#     NumericVector loss(n_splits);
#
#     for (int ii = 0; ii < n_splits; ii++) {
#         double split_value = splits[ii];
#
#         // Vectors to hold left and right indices
#         std::vector<int> left_indices, right_indices;
#
#         // Partition the data based on the split value
#         for (int i = 0; i < n; i++) {
#             if (split_var_values[i] < split_value) {
#                 left_indices.push_back(i);
#             } else {
#                 right_indices.push_back(i);
#             }
#         }
#
#         // If either partition is empty, skip this split
#         if (left_indices.size() == 0 || right_indices.size() == 0) {
#             loss[ii] = NA_REAL;
#             continue;
#         }
#
#         // Extract left and right subsets of X and response
#         arma::mat X_left(left_indices.size(), X.n_cols);
#         arma::vec response_left(left_indices.size());
#         arma::mat X_right(right_indices.size(), X.n_cols);
#         arma::vec response_right(right_indices.size());
#
#         for (size_t i = 0; i < left_indices.size(); i++) {
#             X_left.row(i) = X.row(left_indices[i]);
#             response_left[i] = response[left_indices[i]];
#         }
#         for (size_t i = 0; i < right_indices.size(); i++) {
#             X_right.row(i) = X.row(right_indices[i]);
#             response_right[i] = response[right_indices[i]];
#         }
#
#         // Fit the models using QR decomposition (equivalent to .lm.fit in R)
#         arma::vec coef_left = arma::solve(X_left, response_left);
#         arma::vec coef_right = arma::solve(X_right, response_right);
#
#         // Calculate residuals
#         arma::vec residuals_left = response_left - X_left * coef_left;
#         arma::vec residuals_right = response_right - X_right * coef_right;
#
#         // Compute the sum of squared residuals
#         double rss_left = arma::dot(residuals_left, residuals_left);
#         double rss_right = arma::dot(residuals_right, residuals_right);
#
#         // Store the sum of residuals squared for the split
#         loss[ii] = rss_left + rss_right;
#     }
#
#     return loss;
# }')


splitting_criterion_honest <- function(split.var, formula) {

  split_var_values <- parent.frame()$this_data_fit[[split.var]]
  splits <- sort(unique(split_var_values))
  splits <- unique(quantile(splits, seq(0.05, 0.95, 0.01), type = 1))

  model_frame <- model.frame(formula, parent.frame()$this_data_fit) # for faster .lm.fit
  response <- model_frame[, 1]
  X <- cbind(as.matrix(model_frame[, -1]))

  # Call the Rcpp function
  loss <- splitting_criterion_honest_rcpp(split_var_values, X, response, splits)

  split_at <- splits[which.min(loss)]

  split_var_values <- parent.frame()$this_data_honest[[split.var]]
  model_frame <- model.frame(formula, parent.frame()$this_data_honest) # for faster .lm.fit
  response <- model_frame[, 1]
  X <- cbind(as.matrix(model_frame[, -1]))

  lm_1 <- .lm.fit(X[split_var_values < split_at, , drop = FALSE], response[split_var_values < split_at])$coefficients
  lm_2 <- .lm.fit(X[split_var_values >= split_at, , drop = FALSE], response[split_var_values >= split_at])$coefficients

  return(list(sse = min(loss, na.rm = TRUE), split = split_at, lm_1 = lm_1, lm_2 = lm_2))
}


# # This is the splitting criterion we minimize (SSE [Sum Of Squared Errors]):
# # $SSE = \sum_{i \in S_1} (y_i - \bar(y)1)^2 + \sum_{i \in S_2} (y_i - \bar(y)2)^2$
splitting_criterion_honest_old <- function(split.var,
                                       # data.fit, data.honest,
                                       formula) {

  # browser()

  split_var_values <- parent.frame()$this_data_fit[[split.var]]
  splits <- sort(unique(split_var_values))
  splits <- unique(quantile(splits, seq(0.05, 0.95, 0.01), type = 1))

  model_frame<- model.frame(formula, parent.frame()$this_data_fit) # for faster .lm.fit
  response <- model_frame[,1]
  X <- cbind(as.matrix(model_frame[,-1]))

  loss <- rep(NA, times = length(splits))

  index <- as.numeric(cut(split_var_values, c(-Inf, splits, Inf)))

  for (ii in 1:length(splits)) {

    index_run <- (index <= ii)
    a <- index_run
    b <- !index_run

    lm_1 <- .lm.fit(X[a, , drop = FALSE], response[a])$residuals
    lm_2 <- .lm.fit(X[b, , drop = FALSE], response[b])$residuals

    loss[ii] <- sum(lm_1^2) + sum(lm_2^2)

    # sum((lm_1 - data[split < sp, ]$rv)^2) + sum((lm_2 - data[split >= sp, ]$rv)^2)
    # sum(qlike(pmax(lm_1, min(data$rv)), data[split < sp, ]$rv)) + sum(qlike(pmax(lm_2, min(data$rv)), data[split >= sp, ]$rv))
  }
  # browser()
  # plot(loss)
  split_at <- splits[which.min(loss)]

  split_var_values <- parent.frame()$this_data_honest[[split.var]]
  model_frame<- model.frame(formula, parent.frame()$this_data_honest) # for faster .lm.fit
  response <- model_frame[,1]
  X <- cbind(as.matrix(model_frame[,-1]))

  lm_1 <- .lm.fit(X[split_var_values  < split_at, , drop = FALSE], response[split_var_values < split_at])$coefficients
  lm_2 <- .lm.fit(X[split_var_values  >= split_at, , drop = FALSE], response[split_var_values >= split_at])$coefficients

  return(list(sse = min(loss, na.rm = TRUE), split = split_at, lm_1 = lm_1, lm_2 = lm_2))
}

#' @export
reg_tree_honest <- function(data, formula, split.vars, minsize, mtry = 1/3, data.predict, cl = NULL) {

  # coerce to data.frame
  data <- as.data.frame(data)

  sampled_dates_honest <- sample(unique(data$date), round(length(unique(data$date)) / 2), replace = FALSE)

  data_fit <- data.table(data[data$date %in% sampled_dates_honest, ])
  data_honest <- data.table(data[!(data$date %in% sampled_dates_honest), ])

  # X <- as.matrix(data_fit[, split.vars, drop = FALSE])

  # initialize while loop
  do_splits <- TRUE

  # create output data.frame with splitting rules and observations
  # tree_info <- data.frame(NODE = 1, NOBS = nrow(data), FILTER = NA, TERMINAL = "SPLIT",
  #                         beta0 = NA, beta1 = NA, beta2 = NA, beta3 = NA,
  #                         stringsAsFactors = FALSE)
  tree_info <- data.frame(NODE = 1, NOBS = nrow(data), FILTER = NA, TERMINAL = "SPLIT",
                          beta1 = NA, beta2 = NA, beta3 = NA,
                          stringsAsFactors = FALSE)

  # keep splitting until there are only leafs left
  while (do_splits) {

    # which parents have to be splitted
    to_calculate <- which(tree_info$TERMINAL == "SPLIT")

    for (j in to_calculate) {

      # handle root node
      if (!is.na(tree_info[j, "FILTER"])) {
        # subset data according to the filter
        this_data_fit <- data_fit[eval(parse(text = tree_info[j, "FILTER"]))]
        this_data_honest <- data_honest[eval(parse(text = tree_info[j, "FILTER"]))]
        # get the design matrix
        # X <- as.matrix(this_data_fit[, split.vars, drop = FALSE])

        # if (dim(X) < 120) {
        #   browser()
        # }
      } else {
        this_data_fit <- data_fit
        this_data_honest <- data_honest
      }

      sample_split_vars <- sample(split.vars, mtry)

      splitting <- foreach (kk = sample_split_vars) %do% {
        # browser()
        splitting_criterion_honest(kk,
                                   # data.fit = this_data_fit,
                                   # data.honest = this_data_honest,
                                   formula = formula)
      }

      names(splitting) <- sample_split_vars

      #
      # get the min SSE

      # as.numeric and names workaround only needed for one sorting variable because of
      # different sub-setting behavior of one-dimensional vs. higher dimensional objects
      tmp_splitter <- which.min(map_dbl(splitting, 1))
      names(tmp_splitter) <- names(map_dbl(splitting, 1))[which.min(map_dbl(splitting, 1))]

      # define maxnode
      mn <- max(tree_info$NODE)

      # paste filter rules
      tmp_filter <- c(paste(names(tmp_splitter), ">=",
                            splitting[[tmp_splitter]]$split),
                      paste(names(tmp_splitter), "<",
                            splitting[[tmp_splitter]]$split))

      # Error handling! check if the splitting rule has already been invoked
      split_here  <- !sapply(tmp_filter,
                             FUN = function(x,y) any(grepl(x, x = y)),
                             y = tree_info$FILTER)
      # split_here <- !vapply(tmp_filter,
      #                       FUN = function(pattern) any(grepl(pattern, tree_info$FILTER, perl = TRUE)),
      #                       FUN.VALUE = logical(1))

      # append the splitting rules
      if (!is.na(tree_info[j, "FILTER"])) {
        tmp_filter  <- paste(tree_info[j, "FILTER"],
                             tmp_filter, sep = " & ")
      }

      # get the number of observations in current node
      tmp_nobs <- sapply(tmp_filter,
                         FUN = function(i, x) {
                           nrow(subset(x = x, subset = eval(parse(text = i))))
                         },
                         x = this_data_honest)

      tmp_nobs_fit <- sapply(tmp_filter,
                                FUN = function(i, x) {
                                  nrow(subset(x = x, subset = eval(parse(text = i))))
                                },
                                x = this_data_fit)

      # insufficient minsize for split
      if (any(pmin(c(tmp_nobs
                     , tmp_nobs_fit)) < (minsize))) {
        # browser()
        split_here <- rep(FALSE, 2)
      }

      children_coefficients <- rbind(splitting[[tmp_splitter]]$lm_1, splitting[[tmp_splitter]]$lm_2)
      # colnames(children_coefficients) <- c("beta0", "beta1", "beta2", "beta3")
      colnames(children_coefficients) <- c("beta1", "beta2", "beta3")

      # create children data frame
      children <- data.frame(NODE = c(mn+1, mn+2),
                             NOBS = tmp_nobs,
                             FILTER = tmp_filter,
                             TERMINAL = rep("SPLIT", 2),
                             row.names = NULL)[split_here,]

      if (!all(!split_here)) {
        children <- cbind(children, children_coefficients)
      }

      # overwrite state of current node
      tree_info[j, "TERMINAL"] <- ifelse(all(!split_here), "LEAF", "PARENT")

      # bind everything
      tree_info <- rbind(tree_info, children)

      # check if there are any open splits left
      do_splits <- !all(tree_info$TERMINAL != "SPLIT")
    } # end for
  } # end while

  # browser()


  # calculate fitted values and out-of-sample prediction
  leafs <- tree_info[tree_info$TERMINAL == "LEAF", ]
  predictions <- rep(NA, times = nrow(data.predict))

  if (nrow(leafs) == 1) {
    # stop("No splitting occured. No tree constructed.")
    predictions <- predict(lm(formula, data = data), newdata = data.predict)
  } else {
    data.predict <- as.data.frame(data.predict)
    for (i in seq_len(nrow(leafs))) {
      # ind_prediction <- as.numeric(rownames(subset(as.data.frame(data.predict), eval(parse(text = leafs[i, "FILTER"])))))
      ind_prediction <- attr(subset(data.predict, eval(parse(text = leafs[i, "FILTER"]))), "row.names")
      if (length(ind_prediction) > 0) {
        # prediction_leaf <- unlist(leafs[i, c("beta0", "beta1", "beta2", "beta3")]) %*% t(cbind(1, as.matrix(data.predict[ind_prediction, c("rv_lag_1", "rv_lag_5", "rv_lag_22")])))
        prediction_leaf <- unlist(leafs[i, c("beta1", "beta2", "beta3")]) %*% t(cbind(as.matrix(data.predict[ind_prediction, c("rv_lag_1", "rv_lag_5", "rv_lag_22")])))
        predictions[ind_prediction] <- t(prediction_leaf)
      }
    }
  }

  predictions <-
    data.predict %>%
    select(date, permno, mean_rv) %>%
    mutate(forecast = predictions + mean_rv)

  # return everything
  return(list(tree = tree_info, split_vars = split.vars, predictions = predictions, formula = formula))
}


#' @importFrom purrr map_dbl
#' @importFrom foreach foreach
#' @importFrom foreach '%do%'

#' @export
estimate_har_tree <- function(data, formula, split.vars, minsize, mtry = 2, data.predict) {

  # mtry <- round(length(split.vars) * mtry) # convert to number of subsampled variables

  # coerce to data.frame
  data <- as.data.frame(data)

  sampled_dates_honest <- sample(unique(data$date), round(length(unique(data$date)) / 2), replace = FALSE)

  data_fit <- data[data$date %in% sampled_dates_honest, ]
  data_honest <- data[!(data$date %in% sampled_dates_honest), ]

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
        this_data_fit <- subset(data_fit, eval(parse(text = tree_info[j, "FILTER"])))
        this_data_honest <- subset(data_honest, eval(parse(text = tree_info[j, "FILTER"])))

      } else {
        this_data_fit <- data_fit
        this_data_honest <- data_honest
      }

      sample_split_vars <- sample(split.vars, mtry)

      kk <- NULL
      splitting <- foreach (kk = sample_split_vars) %do% {
        splitting_criterion_honest_package(kk,
                                   data.fit = this_data_fit,
                                   data.honest = this_data_honest,
                                   formula = formula)
      }

      names(splitting) <- sample_split_vars

      # splitting <- list(NA, times = length(sample_split_vars))
      #
      # splitting <- for (kk in sample_split_vars) {
      #   splitting[kk] <- splitting_criterion_honest(kk,
      #                              data.fit = this_data_fit,
      #                              data.honest = this_data_honest,
      #                              formula = formula)
      # }
      #
      # names(splitting) <- sample_split_vars

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

      # insufficient minsize for split
      if (any(tmp_nobs < (minsize))) {
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

  # calculate fitted values and out-of-sample prediction
  leafs <- tree_info[tree_info$TERMINAL == "LEAF", ]
  predictions <- rep(NA, times = nrow(data.predict))

  if (nrow(leafs) == 1) {
    # stop("No splitting occured. No tree constructed.")
    predictions <- predict.lm(lm(formula, data = data), newdata = data.predict)
  } else {
    data.predict <- as.data.frame(data.predict)
    for (i in seq_len(nrow(leafs))) {
      ind_prediction <- attr(subset(data.predict, eval(parse(text = leafs[i, "FILTER"]))), "row.names")
      if (length(ind_prediction) > 0) {
        prediction_leaf <- unlist(leafs[i, c("beta1", "beta2", "beta3")]) %*%
          t(cbind(as.matrix(data.predict[ind_prediction, c("rv_lag_1", "rv_lag_5", "rv_lag_22")])))
        predictions[ind_prediction] <- t(prediction_leaf)
      }
    }
  }
#
#   predictions <-
#     mutate(select(data.predict, date, permno, mean_rv), forecast = predictions + mean_rv)

  df_predictions <-
    data.predict[, c("date", "permno", "mean_rv")]
  df_predictions$forecast <- predictions + df_predictions$mean_rv

  # return everything
  return(list(tree = tree_info, split_vars = split.vars, predictions = df_predictions, formula = formula))
}

# formula, split.vars, minsize, mtry = 1/3, data.predict
#
#
# df_estimation <-
#   rv_panel_data %>%
#   filter(date <= "2004-11-01") %>% # the data is already preaggregated
#   group_by(permno) %>%             # necessary to limit until November
#   mutate(mean_rv = mean(rv_lag_1)) # to avoid look-ahead bias
#
# df_evaluation <-
#   rv_panel_data %>%
#   filter(date >= "2005-01-01") %>%
#   left_join(select(df_estimation, permno, mean_rv) %>% distinct())
#
# estimate_har_tree(
#   df_estimation ,
#   formula = rv_lead_22 ~ 0 + rv_lag_1 + rv_lag_5 + rv_lag_22,
#   split.vars = c("rv_lag_1", "rv_lag_5", "rv_lag_22", "vix_lag"),
#     minsize = 100,
#   mtry = 1/3, # (default)
#   data.predict = df_evaluation
# )
#



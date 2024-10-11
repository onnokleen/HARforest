# Define the Rcpp function
library(Rcpp)
library(RcppArmadillo)

cppFunction('
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
}')

# This is the splitting criterion we minimize (SSE [Sum Of Squared Errors]):
# $SSE = \sum_{i \in S_1} (y_i - \bar(y)1)^2 + \sum_{i \in S_2} (y_i - \bar(y)2)^2$
splitting_criterion_honest_intercept_only <- function(split.var, formula) {
  
  # browser()
  
  # model_frame<- model.frame(formula, data.fit) # for faster .lm.fit
  response <- parent.frame()$this_data_fit[["rv_forecast"]]
  # X <- cbind(as.matrix(model_frame[,-1]))
  
  split_var_values <- parent.frame()$this_data_fit[[split.var]]
  
  # splits <- sort(unique(split_var_values))
  # splits <- unique(quantile(splits, seq(0.05, 0.95, 0.01), type = 1))
  # splits <- sort(unique(quantile(split_var_values, seq(0.05, 0.95, 0.01), type = 1)))
  splits <- unique(quantile(split_var_values, seq(0.05, 0.95, 0.01), type = 1))
  
  # loss <- numeric(length(splits))
  # # browser()
  # 
  # # index <- as.numeric(cut(split_var_values, c(-Inf, splits, Inf)))
  # 
  # for (ii in 1:length(splits)) {
  #   
  #   split_value <- splits[ii]
  #   
  #   index_left <- split_var_values < split_value
  #   index_right <- !index_left
  #   
  #   mean_left <- mean(response[index_left])
  #   mean_right <- mean(response[index_right])
  #   
  #   loss[ii] <- sum((response[index_left] - mean_left)^2) + sum((response[index_right] - mean_right)^2)
  #   
  #   # a <- response[index_left]
  #   # b <- response[index_right]
  #   # 
  #   # loss[ii] <- sum((a - mean(a))^2) + sum((b - mean(b))^2)
  #   
  #   
  #   # sum((lm_1 - data[split < sp, ]$rv)^2) + sum((lm_2 - data[split >= sp, ]$rv)^2)
  #   # sum(qlike(pmax(lm_1, min(data$rv)), data[split < sp, ]$rv)) + sum(qlike(pmax(lm_2, min(data$rv)), data[split >= sp, ]$rv))
  # }
  # # browser()
  # # plot(loss)
  # split_at <- splits[which.min(loss)]
  result <- splitting_criterion_honest_intercept_only_cpp(split_var_values, response, splits)
  
  split_at <- result$split
  
  split_var_values <- parent.frame()$this_data_honest[[split.var]]
  # model_frame<- model.frame(formula, data.honest) # for faster .lm.fit
  response <- parent.frame()$this_data_honest[[toString(formula[[2]])]]
  # X <- cbind(as.matrix(model_frame[,-1]),1)
  
  index <- split_var_values < split_at
  
  lm_1 <- mean(response[index])
  lm_2 <- mean(response[!index])
  
  return(list(sse = min(result$sse, na.rm = TRUE), split = split_at, lm_1 = lm_1, lm_2 = lm_2))
  # return(list(sse = min(loss, na.rm = TRUE), split = split_at, lm_1 = lm_1, lm_2 = lm_2))
}

# # Function to apply complex filtering based on string condition with multiple clauses
# apply_condition <- function(filter_string, data) {
#   # Split the condition string on the '&' symbol to handle multiple conditions
#   conditions <- strsplit(filter_string, " & ")[[1]]
#   
#   # Start with a logical vector of TRUE values, to successively apply conditions
#   logical_result <- rep(TRUE, nrow(data))
#   
#   # Loop over each condition and apply it
#   for (cond in conditions) {
#     # Split each condition into variable, operator, and value
#     parts <- strsplit(cond, " ")[[1]]
#     variable <- parts[1]
#     operator <- parts[2]
#     threshold <- as.numeric(parts[3])
#     
#     # Apply the condition based on the operator
#     if (operator == ">=") {
#       logical_result <- logical_result & (data[[variable]] >= threshold)
#     } else if (operator == "<") {
#       logical_result <- logical_result & (data[[variable]] < threshold)
#     } else {
#       stop("Unsupported operator")
#     }
#   }
#   
#   # Return the number of rows that satisfy all conditions
#   return(nrow(data[logical_result, , drop = FALSE]))
# }
# 
# apply_condition_data <- function(filter_string, data) {
#   conditions <- strsplit(filter_string, " & ")[[1]]
#   logical_result <- rep(TRUE, nrow(data))
#   
#   for (cond in conditions) {
#     parts <- strsplit(cond, " ")[[1]]
#     variable <- parts[1]
#     operator <- parts[2]
#     threshold <- as.numeric(parts[3])
#     
#     if (operator == ">=") {
#       logical_result <- logical_result & (data[[variable]] >= threshold)
#     } else if (operator == "<") {
#       logical_result <- logical_result & (data[[variable]] < threshold)
#     } else {
#       stop("Unsupported operator")
#     }
#   }
#   
#   return(data[logical_result, , drop = FALSE])
# }

reg_tree_honest_intercept_only <- function(data, formula, split.vars, minsize, mtry = 1/3, data.predict, data.oob, cl = NULL) {
  
  # coerce to data.frame
  data <- as.data.frame(data)
  
  sampled_dates_honest <- sample(unique(data$date), round(length(unique(data$date)) / 2), replace = FALSE)
  
  data_fit <- as.data.table(data[data$date %in% sampled_dates_honest, ])
  data_honest <- as.data.table(data[!(data$date %in% sampled_dates_honest), ])
  
  # X <- as.matrix(data_fit[, split.vars, drop = FALSE])
  
  # initialize while loop
  do_splits <- TRUE
  
  # create output data.frame with splitting rules and observations
  # tree_info <- data.frame(NODE = 1, NOBS = nrow(data), FILTER = NA, TERMINAL = "SPLIT", 
  #                         beta0 = NA, beta1 = NA, beta2 = NA, beta3 = NA,
  #                         stringsAsFactors = FALSE)
  tree_info <- data.frame(NODE = 1, NOBS = nrow(data), FILTER = NA, TERMINAL = "SPLIT", 
                          beta0 = NA,
                          stringsAsFactors = FALSE)
  
  # keep splitting until there are only leafs left
  while (do_splits) {
    
    # which parents have to be splitted
    to_calculate <- which(tree_info$TERMINAL == "SPLIT")
    
    for (j in to_calculate) {
      
      # handle root node
      if (!is.na(tree_info[j, "FILTER"])) {
        # subset data according to the filter
        # browser()
        # old
        this_data_fit <- data_fit[eval(parse(text = tree_info[j, "FILTER"]))]
        this_data_honest <- data_honest[eval(parse(text = tree_info[j, "FILTER"]))]
        
        # browser()
        # this_data_fit <- apply_condition_data(tree_info[j, "FILTER"], data_fit)
        # this_data_honest <- apply_condition_data(tree_info[j, "FILTER"], data_honest)
        
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
      
      # browser()
      # splitting <- foreach (kk = sample_split_vars) %do% {
      #   splitting_criterion_honest_intercept_only(kk, 
      #                              # data.fit = this_data_fit,
      #                              # data.honest = this_data_honest,
      #                              formula = formula)
      # }
      # 
      
      # Initialize an empty list to store the results
      splitting <- vector("list", length(sample_split_vars))
      
      # Simple for loop
      for (i in seq_along(sample_split_vars)) {
        kk <- sample_split_vars[i]
        splitting[[i]] <- splitting_criterion_honest_intercept_only(kk, formula = formula)
      }
      
      names(splitting) <- sample_split_vars
      
      list(rep(NULL, times = length(sample_split_vars)))
      
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
      # split_here  <- !sapply(tmp_filter, # old slow code
      #                        FUN = function(x,y) any(grepl(x, x = y)),
      #                        y = tree_info$FILTER)
      
      split_here <- !vapply(tmp_filter, 
                        FUN = function(pattern) any(grepl(pattern, tree_info$FILTER, perl = TRUE)),
                        FUN.VALUE = logical(1))
      
      # append the splitting rules
      if (!is.na(tree_info[j, "FILTER"])) {
        tmp_filter  <- paste(tree_info[j, "FILTER"], 
                             tmp_filter, sep = " & ")
      } 
      # browser()
      
      # get the number of observations in current node
      # old 
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
      if (any(c(tmp_nobs, tmp_nobs_fit) < (minsize))) {
        # browser()
        split_here <- rep(FALSE, 2)
      }
      
      children_coefficients <- rbind(splitting[[tmp_splitter]]$lm_1, splitting[[tmp_splitter]]$lm_2)
      # colnames(children_coefficients) <- c("beta0", "beta1", "beta2", "beta3")
      colnames(children_coefficients) <- c("beta0")
      
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
      # tree_info[j, "TERMINAL"] <- ifelse(all(!split_here), "LEAF", "PARENT") #old
      # Update TERMINAL state for the current node
      if (all(!split_here)) {
        tree_info[j, "TERMINAL"] <- "LEAF"
      } else {
        tree_info[j, "TERMINAL"] <- "PARENT"
      }
      
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
    for (i in seq_len(nrow(leafs))) {
      ind_prediction <- as.numeric(rownames(subset(as.data.frame(data.predict), eval(parse(text = leafs[i, "FILTER"])))))
      if (length(ind_prediction) > 0) {
        # prediction_leaf <- unlist(leafs[i, c("beta0", "beta1", "beta2", "beta3")]) %*% t(cbind(1, as.matrix(data.predict[ind_prediction, c("rv_lag_1", "rv_lag_5", "rv_lag_22")])))
        # prediction_leaf <- unlist(leafs[i, c("beta0")]) %*% t(cbind(as.matrix(data.predict[ind_prediction, c("rv_lag_1", "rv_lag_5", "rv_lag_22")])))
        predictions[ind_prediction] <- leafs[i, c("beta0")]
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

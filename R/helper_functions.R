#' This is the splitting criterion we minimize (i.e., sum of squares)
#'
#' @importFrom stats .lm.fit
#' @importFrom stats lm
#' @importFrom stats model.frame
#' @importFrom stats predict.lm
#' @importFrom stats quantile
#'
#' @keywords internal
splitting_criterion_honest <- function(split.var, data.fit, data.honest, formula) {


  model_frame<- model.frame(formula, data.fit) # for faster .lm.fit
  response <- model_frame[,1]
  X <- cbind(as.matrix(model_frame[,-1]))

  splits <- sort(unique(data.fit[, split.var]))
  splits <- unique(quantile(splits, seq(0.05, 0.95, 0.01), type = 1))

  split_var_values <- data.fit[, split.var]

  model_frame<- model.frame(formula, data.fit) # for faster .lm.fit
  response <- model_frame[,1]
  X <- cbind(as.matrix(model_frame[,-1]))

  loss <- rep(NA, times = length(splits))

  for (ii in 1:length(splits)) {
    sp <- splits[ii]

    lm_1 <- .lm.fit(X[split_var_values < sp, , drop = FALSE], response[split_var_values < sp])$residuals
    lm_2 <- .lm.fit(X[split_var_values >= sp, , drop = FALSE], response[split_var_values >= sp])$residuals

    loss[ii] <- sum(lm_1^2) + sum(lm_2^2)
  }
  split_at <- splits[which.min(loss)]

  lm_1 <- .lm.fit(X[split_var_values  < split_at, , drop = FALSE], response[split_var_values < split_at])$coefficients
  lm_2 <- .lm.fit(X[split_var_values  >= split_at, , drop = FALSE], response[split_var_values  >= split_at])$coefficients

  split_var_values <- data.honest[, split.var]
  model_frame<- model.frame(formula, data.honest) # for faster .lm.fit
  response <- model_frame[,1]
  X <- cbind(as.matrix(model_frame[,-1]))

  lm_1 <- .lm.fit(X[split_var_values  < split_at, , drop = FALSE], response[split_var_values < split_at])$coefficients
  lm_2 <- .lm.fit(X[split_var_values  >= split_at, , drop = FALSE], response[split_var_values >= split_at])$coefficients

  return(list(sse = min(loss, na.rm = TRUE), split = split_at, lm_1 = lm_1, lm_2 = lm_2))
}

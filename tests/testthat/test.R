# library(HARforest)
# library(foreach)
# library(dplyr)



test_that("trigonometric functions match identities", {
  library(dplyr)
  split_vars <- c("rv_lag_1", "rv_lag_5", "rv_lag_22", "vix_lag")

  # Currently, the data has to be structed and named as in our sample data "rv_panel_data"
  # which is part of our package.
  df_estimation <-
    rv_panel_data %>%
    filter(date <= "2004-11-01") %>%     # the data is already preaggregated
    group_by(permno) %>%                 # necessary to limit until November
    mutate(mean_rv = mean(rv_lag_1)) %>% # to avoid look-ahead bias
    mutate(across(c(rv_lead_22, rv_lag_1, rv_lag_5, rv_lag_22), ~ . - mean_rv))

  df_evaluation <-
    rv_panel_data %>%
    filter(date >= "2005-01-01") %>%
    left_join(select(df_estimation, permno, mean_rv) %>% distinct()) %>%
    mutate(across(c(rv_lead_22, rv_lag_1, rv_lag_5, rv_lag_22), ~ . - mean_rv))

  set.seed(123)

  test_tree <- estimate_har_tree(
    df_estimation ,
    formula = rv_lead_22 ~ 0 + rv_lag_1 + rv_lag_5 + rv_lag_22,
    split.vars = split_vars,
    minsize = 100,
    mtry = 1, # (default)
    data.predict = df_evaluation
  )

  expect_equal(test_tree$tree[4, ]$NOBS, 1167)
})

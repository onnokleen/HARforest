# HARforest

<!-- badges: start -->
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![CRAN status](https://www.r-pkg.org/badges/version/HARforest)](https://CRAN.R-project.org/package=HARforest)
<!-- badges: end -->

The goal of HARforest is to enable practitioners and researchers alike to employ the HAR forest developed by Onno Kleen and Anastasija Tetereva (2022, [doi:10.2139/ssrn.4161957](https://10.2139/ssrn.4161957)).

## Installation

You can install the development version of HARforest like this:

``` r
library(devtools)
install_github("onnokleen/HARforest")
```

We are working on making the code available on CRAN soon.

## Example

This is a basic example which shows you how to estimate a HAR tree, how to grow a forest, and how to predict with it:

```{r example}
library(HARforest)
library(foreach)
library(dplyr)

split_vars <- c("rv_lag_1", "rv_lag_5", "rv_lag_22", "vix_lag")

# Currently, the data still has to be structed and named as in our sample data "rv_panel_data"
# We will update the package to be more flexible before releasing it on CRAN
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

# Estimate single HAR tree
estimate_har_tree(
  df_estimation ,
  formula = rv_lead_22 ~ 0 + rv_lag_1 + rv_lag_5 + rv_lag_22,
  split.vars = split_vars,
  minsize = 100,
  mtry = 1/3, # (default)
  data.predict = df_evaluation
)

# Grow multiple trees to form a forest
# Parallelization can be done with doParallel package and using %dopar% instead of %do%

n_trees <- 10

tree_list <- foreach (ii = 1:n_trees) %do% {
  tree <- estimate_har_tree(df_estimation,
                            formula = rv_lead_22 ~ 0 + rv_lag_1 + rv_lag_5 + rv_lag_22 ,
                            split.vars = split_vars,
                            minsize = 100,
                            data.predict = df_evaluation)
  environment(tree$formula) <- NULL
  tree
}

# Aggregate predictions of trees
tree_predictions <-
  tree_list %>%
  lapply(., function(x) x$predictions$forecast) %>%
  do.call(cbind, .) %>%
  rowMeans()

# Align them with evaluation sample
df_evaluation %>%
  mutate(har_forest = tree_predictions)
  
```


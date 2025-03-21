# HARforest

## Overview

HARforest is an R package that implements the heterogeneous autoregressive (HAR) forest methodology for volatility forecasting. This innovative approach combines the simplicity of HAR models with the power of random forests to capture complex, non-linear relationships in financial time series.

### Key Features

- **Panel Data Support**: Handles multiple time series simultaneously, making it ideal for portfolio analysis
- **Honest Estimation**: Implements honest estimation to prevent overfitting
- **Easy Integration**: Seamlessly works with standard R data structures and the tidyverse ecosystem

### Use Cases

- Volatility forecasting for financial assets
- Risk management and portfolio optimization
- Market regime detection
- High-frequency trading applications
- Academic research in financial econometrics

The methodology is based on the paper "A forest full of risk forecasts for managing volatility" by Kleen, Santos, and Tetereva (2025, [doi:10.2139/ssrn.4161957](https://10.2139/ssrn.4161957)).

## Installation

You can install the development version of HARforest like this:

``` r
library(devtools)
install_github("onnokleen/HARforest")
```

## Quick Start Guide

### Prerequisites

Before using HARforest, make sure you have the following packages installed:
``` r
install.packages(c("dplyr", "foreach", "purrr"))
```

### Basic Usage

Here's a step-by-step example showing how to:
1. Prepare your data
2. Estimate a single HAR tree
3. Grow a forest of trees
4. Make predictions

``` r
# Load required packages
library(HARforest)
library(foreach)
library(dplyr)

# Define splitting variables for the tree
split_vars <- c("rv_lag_1", "rv_lag_5", "rv_lag_22", "vix_lag")

# Prepare estimation data
# Note: Your data should follow the same structure as the sample data "rv_panel_data"
df_estimation <-
  rv_panel_data %>%
  filter(date <= "2004-11-01") %>%     # Split data into estimation period
  group_by(permno) %>%                 # Group by asset identifier
  mutate(mean_rv = mean(rv_lag_1)) %>% # Calculate mean realized volatility
  mutate(across(c(rv_lead_22, rv_lag_1, rv_lag_5, rv_lag_22), ~ . - mean_rv)) # Center variables

# Prepare evaluation data
df_evaluation <-
  rv_panel_data %>%
  filter(date >= "2005-01-01") %>%     # Split data into evaluation period
  left_join(select(df_estimation, permno, mean_rv) %>% distinct()) %>%
  mutate(across(c(rv_lead_22, rv_lag_1, rv_lag_5, rv_lag_22), ~ . - mean_rv))

# Step 1: Estimate a single HAR tree
tree <- estimate_har_tree(
  df_estimation,
  formula = rv_lead_22 ~ 0 + rv_lag_1 + rv_lag_5 + rv_lag_22,
  split.vars = split_vars,
  minsize = 100,    # Minimum number of observations in each leaf
  mtry = 1/3,       # Number of variables to consider at each split
  data.predict = df_evaluation
)

# Step 2: Grow a forest of trees
# Note: For parallel processing, use doParallel package and %dopar% instead of %do%
n_trees <- 10

tree_list <- foreach(ii = 1:n_trees) %do% {
  tree <- estimate_har_tree(
    df_estimation,
    formula = rv_lead_22 ~ 0 + rv_lag_1 + rv_lag_5 + rv_lag_22,
    split.vars = split_vars,
    minsize = 100,
    data.predict = df_evaluation
  )
  environment(tree$formula) <- NULL
  tree
}

# Step 3: Aggregate predictions from all trees
har_forest_predictions <-
  tree_list %>%
  lapply(., function(x) x$predictions$forecast) %>%
  do.call(cbind, .) %>%
  rowMeans()

# Step 4: Combine predictions with evaluation data
results <- df_evaluation %>%
  mutate(har_forest = har_forest_predictions)
```

### Important Notes

- The data should be structured as a panel dataset with columns for dates, asset identifiers, and realized volatility measures
- Variables should be named consistently with the sample data structure
- Consider using parallel processing for large forests
- The `minsize` parameter controls the minimum number of observations in each leaf node
- The `mtry` parameter determines the number of variables considered at each split



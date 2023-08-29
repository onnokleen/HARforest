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

This is a basic example which shows you how to estimate a HAR forest and to predict with it:

```{r example}
library(HARforest)
df_estimation <-
  rv_panel_data %>%
  filter(date <= "2004-11-01") %>% # the data is already preaggregated
  group_by(permno) %>%             # necessary to limit until November
  mutate(mean_rv = mean(rv_lag_1)) # to avoid look-ahead bias

df_evaluation <-
  rv_panel_data %>%
  filter(date >= "2005-01-01") %>%
  left_join(select(df_estimation, permno, mean_rv) %>% distinct())

estimate_har_tree(
  df_estimation ,
  formula = rv_lead_22 ~ 0 + rv_lag_1 + rv_lag_5 + rv_lag_22,
  split.vars = c("rv_lag_1", "rv_lag_5", "rv_lag_22", "vix_lag"),
  minsize = 100,
  mtry = 1/3, # (default)
  data.predict = df_evaluation
)
```


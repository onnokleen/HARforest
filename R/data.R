#' Stock returns and financial conditions.
#'
#' A dataset containing the S&P 500 stock returns and the NFCI
#'
#' @format A data frame with 16,004 rows and 7 variables:
#' \describe{
#'   \item{date}{date}
#'   \item{permno}{a permanent stock identifier}
#'   \item{rv_lead_22}{5-minute realized variances}
#'   \item{rv_lag_1}{a dummy for each year/week combination}
#'   \item{rv_lag_5}{National Financial Conditions Index}
#'   \item{rv_lag_22}{National Financial Conditions Index}
#'   \item{vix_lag}{National Financial Conditions Index}#' }
"rv_panel_data"

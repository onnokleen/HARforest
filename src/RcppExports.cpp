// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// splitting_criterion_honest_rcpp
Rcpp::NumericVector splitting_criterion_honest_rcpp(const arma::vec& split_var_values, const arma::mat& X, const arma::vec& response, const Rcpp::NumericVector& splits);
RcppExport SEXP _HARforest_splitting_criterion_honest_rcpp(SEXP split_var_valuesSEXP, SEXP XSEXP, SEXP responseSEXP, SEXP splitsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type split_var_values(split_var_valuesSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type response(responseSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type splits(splitsSEXP);
    rcpp_result_gen = Rcpp::wrap(splitting_criterion_honest_rcpp(split_var_values, X, response, splits));
    return rcpp_result_gen;
END_RCPP
}
// splitting_criterion_honest_intercept_only_cpp
List splitting_criterion_honest_intercept_only_cpp(NumericVector split_var_values, NumericVector response, NumericVector splits);
RcppExport SEXP _HARforest_splitting_criterion_honest_intercept_only_cpp(SEXP split_var_valuesSEXP, SEXP responseSEXP, SEXP splitsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type split_var_values(split_var_valuesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type response(responseSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type splits(splitsSEXP);
    rcpp_result_gen = Rcpp::wrap(splitting_criterion_honest_intercept_only_cpp(split_var_values, response, splits));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_HARforest_splitting_criterion_honest_rcpp", (DL_FUNC) &_HARforest_splitting_criterion_honest_rcpp, 4},
    {"_HARforest_splitting_criterion_honest_intercept_only_cpp", (DL_FUNC) &_HARforest_splitting_criterion_honest_intercept_only_cpp, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_HARforest(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
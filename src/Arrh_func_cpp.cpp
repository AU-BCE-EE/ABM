#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector Arrh_func_cpp(NumericVector A, NumericVector E, double R, double temp_K) {
  return A * exp(-E / (R * temp_K));
}
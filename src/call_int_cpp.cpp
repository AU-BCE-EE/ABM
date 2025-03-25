#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector call_int(Rcpp::Function temp_C_fun, double x) {
  // Call the custom R function passed as an argument
  NumericVector temp_C = temp_C_fun(x);
  
  // Return the result
  return temp_C;
}

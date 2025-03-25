#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector extract_xa_cpp(NumericVector y, int n_mic) {
  // Subset y from 1 to n_mic (in C++, indexing starts from 0)
  NumericVector xa = y[Range(0, n_mic - 1)];
  
  return xa;
}
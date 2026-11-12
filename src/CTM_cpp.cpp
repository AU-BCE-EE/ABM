#include <Rcpp.h>
#include "ctm.h"

using namespace Rcpp;

// Cardinal temperature model, vectorized over microbial groups
// Only first element of tt is used
// [[Rcpp::export(rng = false)]]
NumericVector CTM_cpp(NumericVector tt, NumericVector top, NumericVector tmin,
                      NumericVector tmax, NumericVector yopt) {
  int n = top.size();
  NumericVector y(n);

  for (int i = 0; i < n; i++) {
    y[i] = ctm1(tt[0], top[i], tmin[i], tmax[i], yopt[i]);
  }

  return y;
}

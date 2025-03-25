#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector Arrh_func_cpp(NumericVector A, NumericVector E, double R, double temp_K, 
                            double scale_alpha_opt, double alpha_opt_scale_type, double alpha_opt_scale_CP
                            ) {
  
  // index 7 is urea, index 3 and 4 is CPs and CPf
  
  NumericVector y = A * exp(-E / (R * temp_K));
  y[Rcpp::Range(0,6)] = y[Rcpp::Range(0,6)] * scale_alpha_opt * alpha_opt_scale_type;
  y[Rcpp::Range(3,4)] = y[Rcpp::Range(3,4)] * alpha_opt_scale_CP;
  
  return y;
}
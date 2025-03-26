#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector CTM_cpp(NumericVector tt, NumericVector top, NumericVector tmin, 
                            NumericVector tmax, NumericVector yopt) {
  int n = top.size();
  NumericVector y(n);
  
    for (int i = 0; i < n; i++) {
      
      if (top[i] - tmin[i] < (tmax[i] - tmin[i])/2) {
        
          y[i] = yopt[i] * ((tt[0] - tmin[i]) * pow(tt[0] - tmax[i], 2)) / 
          ((top[i] - tmax[i]) * ((top[i] - tmax[i]) * (tt[0] - top[i]) - 
          (top[i] - tmin[i]) * (top[i] + tmax[i] - 2*tt[0])));
        
        if (y[i] < 0){
          y[i] = 0;
        }
        if (tt[0] <= tmin[i] || tt[0] >= tmax[i]) {
          y[i] = 0;
        }
        
      } else{
         
          y[i] = yopt[i] * ((tt[0] - tmax[i]) * pow(tt[0] - tmin[i], 2)) / 
          ((top[i] - tmin[i]) * ((top[i] - tmin[i]) * (tt[0] - top[i]) - 
          (top[i] - tmax[i]) * (top[i] + tmin[i] - 2*tt[0])));
        
        if (y[i] < 0){
          y[i] = 0;
        }
        if (tt[0] <= tmin[i] || tt[0] >= tmax[i]) {
          y[i] = 0;
        }
      }
    }

  return y; 

}
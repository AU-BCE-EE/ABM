#include <Rcpp.h>
#include <string>

using namespace Rcpp;

// [[Rcpp::export]]
LogicalVector grepl_cpp(std::string pattern, CharacterVector x) {
  LogicalVector result(x.size());
  
  for (R_xlen_t i = 0; i < x.size(); i++) {
    std::string str = Rcpp::as<std::string>(x[i]);
    
    if (x[i] == NA_STRING) {
      result[i] = NA_LOGICAL;
    } else {
      // Use find() to check if the pattern occurs at the start (position 0)
      if (str.find(pattern) == 0) {
        result[i] = true;
      } else {
        result[i] = false;
      }
    }
  }
  
  return result;
}

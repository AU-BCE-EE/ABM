rm(list = ls())


library(Rcpp)
library(microbenchmark)

# generate some temperature data
temp_C <- data.frame(time = 1:365, temp_C = 15 + 5* sin(1:365))

# use current ABM approach where a R approx_fun() is made based on the temp data
source('../R/makeTimeFunc.R')

temp_C_fun <- makeTimeFunc(temp_C, approx_method = 'linear')

# make a C++ function that interpolates at given x values
cppFunction('NumericVector temp_C_fun_cpp(NumericVector t, NumericVector y, NumericVector new_t) {
  
  int n = new_t.size();  // Number of new time points
  int m = t.size();      // Number of existing time points
  NumericVector temp_C(n);  // Result vector for the interpolated values

  // Access data directly using pointers
  const double* t_ptr = t.begin();
  const double* y_ptr = y.begin();
  const double* new_t_ptr = new_t.begin();
  double* temp_C_ptr = temp_C.begin();

  // Binary search approach for each new_t
  for (int i = 0; i < n; i++) {
    double t_query = new_t_ptr[i];
    
    // Use std::lower_bound to find the first element that is not less than t_query
    auto it = std::lower_bound(t_ptr, t_ptr + m, t_query);
    int idx = it - t_ptr;  // Get the index of the found position
    
    // Boundary conditions: if t_query is exactly t[idx], no interpolation needed
    if (idx < m && t_ptr[idx] == t_query) {
      temp_C_ptr[i] = y_ptr[idx];  // Exact match
    } else if (idx == 0) {
      temp_C_ptr[i] = y_ptr[0];  // Extrapolate to the left (before the first value)
    } else if (idx == m) {
      temp_C_ptr[i] = y_ptr[m - 1];  // Extrapolate to the right (after the last value)
    } else {
      // Linear interpolation between t[idx-1] and t[idx]
      double t1 = t_ptr[idx - 1], t2 = t_ptr[idx];
      double y1 = y_ptr[idx - 1], y2 = y_ptr[idx];
      temp_C_ptr[i] = y1 + (t_query - t1) * (y2 - y1) / (t2 - t1);  // Interpolation formula
    }
  }

  return temp_C;
}')

new_t_short <- runif(100, 1, 365)
cpp_out <- temp_C_fun_cpp(t = c(1:365), y = temp_C$temp_C, new_t = new_t_short)

R_out <- temp_C_fun(new_t_short)
R_approx <- function(x, y, new_t){
  return (approx(x, y, new_t)$y)
}

R_approx_out <- R_approx(x = c(1:365), y = temp_C$temp_C, new_t = new_t_short)

# are same results produced?
all(cpp_out == R_out)
all(cpp_out == R_approx_out)

# time of different functions
microbenchmark::microbenchmark(cpp_time = temp_C_fun_cpp(t = c(1:365), y = temp_C$temp_C, new_t = new_t_short),
                               R_time = temp_C_fun(new_t_short),
                               R_approx_time = R_approx(x = c(1:365), y = temp_C$temp_C, new_t = new_t_short))

# with many interpolations how does speed performance change?
new_t_long <- runif(100000, 1, 365)

microbenchmark::microbenchmark(cpp_time = temp_C_fun_cpp(t = c(1:365), y = temp_C$temp_C, new_t = new_t_long),
                               R_time = temp_C_fun(new_t_long),
                               R_approx_time = R_approx(x = c(1:365), y = temp_C$temp_C, new_t = new_t_long))

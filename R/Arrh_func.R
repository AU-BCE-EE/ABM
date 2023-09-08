Arrh_func <- function(A, E, R, temp_K){
  
  y <- A * exp(-E/(R * temp_K))
  
  return(y)
}
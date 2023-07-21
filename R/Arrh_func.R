Arrh_func <- function(lnA, E, R, temp_K){
  
  y <- lnA * exp(-E/(R * temp_K))
  
  return(y)
}
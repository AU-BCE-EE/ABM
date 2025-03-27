rates_test <- function(t, y, parms, temp_C_fun = temp_C_fun, pH_fun = pH_fun, 
                  SO4_inhibition_fun = SO4_inhibition_fun, 
                  conc_fresh_fun = conc_fresh_fun, xa_fresh_fun = xa_fresh_fun){

out <- rates_cpp(t, y, parms, temp_C_fun, pH_fun, SO4_inhibition_fun, conc_fresh_fun, 
          xa_fresh_fun, CTM_cpp)  
  
list2env(out, environment())  

  return(list(derivatives, rates))
  
}  
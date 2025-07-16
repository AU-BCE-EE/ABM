# Sorts out slurry production values and removal timing

calcProdRem <- function(var_pars, approx_method) {
  
  # Determine slurry removal quantity in each time interval
  # Note final 0--alignment is a bit tricky
  if (approx_method %in% c('late', 'mid')) {
    removals <- - c(0, diff(var_pars$var[-nrow(var_pars$var), 'slurry_mass']), 0) > 0
  } else if (approx_method == 'early') {
    removals <- - c(diff(var_pars$var[, 'slurry_mass']), 0) > 0
  } 
  var_pars$var$removal <- removals

  # Add column with slurry production rate
  slurry_prod_rate_t <- c(diff(var_pars$var[, 'slurry_mass']) / diff(var_pars$var[, 'time']), 0) 
  slurry_prod_rate_t[slurry_prod_rate_t < 0] <- 0
  slurry_prod_rate_t[!is.finite(slurry_prod_rate_t)] <- 0
  var_pars$var$slurry_prod_rate <- slurry_prod_rate_t

  return(var_pars)

}

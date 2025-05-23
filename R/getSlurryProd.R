# Extracts slurry production rate over time

getSlurryProd <- function(pars) {
  
  slurry_prod_rate_t <- c(0, diff(pars$slurry_mass[, 'slurry_mass']) / diff(pars$slurry_mass[, 'time'])) 
  slurry_prod_rate_t[slurry_prod_rate_t < 0] <- 0
  slurry_prod_rate_t[!is.finite(slurry_prod_rate_t)] <- 0

  return(slurry_prod_rate_t)
}


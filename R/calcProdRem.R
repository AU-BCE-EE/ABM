# Sorts out slurry production values and removal timing

calcProdRem <- function(pars) {
  
  if (is.null(pars$var)) {
    return(pars)
  }

  # Removals ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Note final 0--alignment is a bit tricky
  if (pars$approx_method %in% c('late', 'mid')) {
    removals <- - c(0, diff(pars$var[-nrow(pars$var), 'slurry_mass']), 0) > 0
  } else if (pars$approx_method == 'early') {
    removals <- - c(diff(pars$var[, 'slurry_mass']), 0) > 0
  } 
  pars$var$removal <- removals

  # Slurry production rate ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  slurry_prod_rate_t <- c(diff(pars$var[, 'slurry_mass']) / diff(pars$var[, 'time']), 0) 
  slurry_prod_rate_t[slurry_prod_rate_t < 0] <- 0
  slurry_prod_rate_t[!is.finite(slurry_prod_rate_t)] <- 0
  pars$var$slurry_prod_rate <- slurry_prod_rate_t

  # Residual slurry for emptying ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if (pars$approx_method == 'late') {
    pars$var$resid_mass <- pars$var$slurry_mass
  } else {
    pars$var$resid_mass <- c(pars$var$slurry_mass[-1], 0)
  }

  return(pars)
}

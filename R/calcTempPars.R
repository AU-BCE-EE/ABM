
calcTempPars <- function(pars, y) {
  
  pars$temp_K <- pars$temp_C + 273.15
  
  # Temperature-dependent derivative vectors
  pars$alpha <- pars$qhat <- 0 * y
  # Hydrolysis rate (vectorized)
  pars$alpha[pars$subs] <- CTM_cpp(pars$temp_K, pars$T_opt_hyd, pars$T_min_hyd, pars$T_max_hyd, pars$hydrol_opt)
  # Microbial substrate utilization rate (vectorized calculation)
  pars$qhat[pars$grps] <- CTM_cpp(pars$temp_K, pars$T_opt, pars$T_min, pars$T_max, pars$qhat_opt)

  return(pars)

}

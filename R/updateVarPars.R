# Copy time-variable parameters into their normal par position for an interval

updateVarPars <- function(pars, y, i) {

  if (any(names(pars$var) != 'time')) {
    vdat <- pars$var[, names(pars$var) != 'time', drop = FALSE]

    for (j in 1:ncol(vdat)) {
      pn <- names(vdat)[j]
      newval <- vdat[i, j] 
      if (is.list(newval)) {
        newval <- newval[[1]]
      }
      pars[[pn]] <- newval
    }

    # Convert temperature to K in case any C values were given in var
    pars <- tempsC2K(pars, cutoff = 200)
    pars$temp_K <- pars$temp_C + 273.15
  }
  
  # Temperature-dependent derivative vectors
  pars$alpha <- pars$qhat <- 0 * y
  # Hydrolysis rate (vectorized)
  pars$alpha[pars$subs] <- CTM_cpp(pars$temp_K, pars$T_opt_hyd, pars$T_min_hyd, pars$T_max_hyd, pars$hydrol_opt)
  # Microbial substrate utilization rate (vectorized calculation)
  pars$qhat[pars$grps] <- CTM_cpp(pars$temp_K, pars$T_opt, pars$T_min, pars$T_max, pars$qhat_opt)

  return(pars)

}

# Creates list of interval times for simulation

makeTimeList <- function(pars, times, days, delta_t) {
  
  # Sort out times returned by ODE solver
  if (is.null(times)) {
    times <- seq(0, days, by = delta_t)
  }

  # Include days argument in times vector
  times <- sort(unique(c(times, pars$var$time)))
    
  # Note about time: 1) All simulations start at 0, 2) days must be at least as long as var data
  # Note that this works even with t_end = NULL (is ignored)
  # Note the "dummy" placeholder in position 1 (and extra + 1 in n_int)
  n_int <- nrow(pars$var)
  st <- pars$var[, 'time']
  timelist <- cumtime <- as.list(rep(0, n_int))
  for (i in 2:n_int) {
    tt <- times[times > st[i - 1] & times <= st[i]] # slow speed of this line
    # Simulation intervals should end exactly at emptying time, so it is added here through st[i] for cases where there is not alignment
    tt <- unique(c(tt, st[i]))
    if (length(tt) == 0) { 
      # Not clear when this might happen, but it would be bad
      stop('No times for interval (row) ', i, ' in time list for some reason. xyn917')
    } else {
      cumtime[[i]] <- max(tt)
      tt <- tt - cumtime[[i - 1]]
      tt <- unique(c(0, tt))
    }
    timelist[[i]] <- tt
  }

  return(timelist)

}

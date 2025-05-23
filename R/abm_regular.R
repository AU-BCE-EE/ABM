abm_regular <- function(
      days, 
			delta_t, 
			times_regular, 
			y, 
			pars, 
			starting = NULL, 
			temp_C_fun = temp_C_fun, 
			pH_fun = pH_fun) { 

  # Figure out timing
  intervals <- getRegTimes(pars, days)

  # Empty data frame for holding results
  dat <- NULL

  # Time trackers
  # Time remaining to run
  t_rem <- days
  # Time that has already run
  t_run <- 0

  # Start the time (emptying) loop
  for (i in 1:intervals$n_int) {

    # Sort out call duration
    t_call <- intervals$t_int[i]
    
    # Need some care with times to make sure t_call is last one in case it is not multiple of delta_t
    times <- sort(unique(round(c(seq(0, t_call, by = min(t_rem, delta_t)), t_call), 5)))
    
    # Add run time to pars so rates() can use actual time to calculate temp_C and pH
    pars$t_run <- t_run
    pars$t_call <- t_call
    pars$times <- times

    # Call up ODE solver
    out <- deSolve::lsoda(y = y, 
                          times = times, 
                          rates, 
                          parms = pars, 
                          temp_C_fun = temp_C_fun, 
                          pH_fun = pH_fun)

    # Extract new state variable vector from last row
    y <- getLastState(out, y)
    
    # Empty channel (instantaneous changes at end of day) in preparation for next lsoda call
    y <- emptyStore(y, resid_mass = pars$resid_mass, resid_enrich = pars$resid_enrich)
    y.eff <- y$eff
    y <- y$store
   
    # Washing, increase slurry mass
    if (intervals$wash[i]) {
      y['slurry_mass'] <- y['slurry_mass'] + pars$wash_water

      # And empty again, leaving same residual mass as always, and enriching for particles
      y <- emptyStore(y, resid_mass = pars$resid_mass, resid_enrich = pars$resid_enrich)
      y.eff <- y$eff
      y <- y$store

      # Run for rest period with no influent 
      if (pars$rest_d > 0) {
        times <- seq(0, pars$rest_d, delta_t)
        parsr <- pars
        parsr$slurry_prod_rate <- 0
        outr <- deSolve::lsoda(y = y, 
                              times = times, 
                              rates, 
                              parms = parsr, 
                              temp_C_fun = temp_C_fun, 
                              pH_fun = pH_fun)
        # Extract new state variable vector from last row
        y <- getLastState(outr, y)
        
        # Correct time in outr and combine with main output
        out <- addOut(main = out, new = outr, t_add = out[nrow(out), 'time'])
      }
    }
    
    # Clean up and stack output with earlier results
    dat <- addOut(main = dat, new = out, t_add = t_run, y.eff = y.eff)
    
    # Update time remaining and time run so far
    t_rem <- t_rem - t_call - intervals$wash[i] * pars$rest_d
    t_run <- t_run + t_call + intervals$wash[i] * pars$rest_d

  }

  if (!is.null(times_regular)) {
    times_regular <- sort(unique(c(times_regular, days)))
    dat <- dat[dat$time %in% times_regular, ]
  }

  return(dat)

}

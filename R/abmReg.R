abmReg <- function(
      days, 
			delta_t, 
			times_regular, 
			y, 
			pars) { 

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

    # Calculate log Ka for speciation (temperature dependent)
    pars <- calcKa(pars)
  
    # Add run time to pars so rates() can use actual time to calculate temp_C and pH
    pars$t_run <- t_run
    pars$t_call <- t_call
    pars$times <- times
    
    # Call up ODE solver
    out <- deSolve::lsoda(y = y, 
                          times = times, 
                          rates, 
                          parms = pars)
    
    # Update time remaining and time run so far
    t_run <- t_run + t_call

    # Extract new state variable vector from last row
    y <- getLastState(out, y)
      
    # See if emptying is needed
    # As long as simulation is not out of time, the end of an interval means emptying occurs
    if (t_run < days) {

      # Empty channel (instantaneous changes at end of day) in preparation for next lsoda call
      y <- emptyStore(y, resid_mass = pars$resid_mass, resid_enrich = pars$resid_enrich)
      y.eff <- y$eff
      y <- y$store

      # Add effluent to output
      out <- addEff(out, y.eff)
   
      # Optional washing
      if (intervals$wash[i]) {
        
        y['slurry_mass'] <- y['slurry_mass'] + pars$wash_water

        # And empty again, leaving same residual mass as always, and enriching for particles
        y <- emptyStore(y, resid_mass = pars$resid_mass, resid_enrich = pars$resid_enrich)
        y.eff <- y$eff
        y <- y$store

        # Add post-wash state to output data frame
        out <- addOut(main = out, new = y)
        
        # Add washing effluent to output
        out <- addEff(out, y.eff)

      }
      
      # Run for rest period with no influent 
      if (pars$rest_d > 0) {
        times <- seq(0, pars$rest_d, delta_t)
        parsr <- pars
        parsr$slurry_prod_rate <- 0
        outr <- deSolve::lsoda(y = y, 
                              times = times, 
                              rates, 
                              parms = parsr)
        # Extract new state variable vector from last row
        y <- getLastState(outr, y)
        
        # Add effluent to output
        outr <- addEff(outr, 0 * y.eff)

        # Correct time in outr and combine with main output
        out <- addOut(main = out, new = outr)
      
        # Add rest time to t_run
        t_run <- t_run + intervals$wash[i] * pars$rest_d
      }
      
    } else {
      
      # Make sure there are some effluent columns (with 0)
      y.eff <- 0 * emptyStore(y)$eff
      out <- addEff(out, y.eff)
  
    }
    
    # Stack output with earlier results
    dat <- addOut(main = dat, new = out)
    
    # Update time
    t_rem <- t_rem - t_call - intervals$wash[i] * pars$rest_d
    
  }

  # Subset to regular times
  if (!is.null(times_regular)) {
    times_regular <- sort(unique(c(times_regular, days)))
    dat <- dat[dat$time %in% times_regular, ]
  }

  return(dat)

}

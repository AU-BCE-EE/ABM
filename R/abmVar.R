abmVar <-
  function(days, 
           delta_t, 
           times, 
           y, 
           pars, 
           warn) {
    
  pars$abm_regular <- FALSE
  
  # Extract washing mass
  wash_water <- getWashWater(pars)
  
  # Extract slurry_mass vector for use in emptying calculations
  slurry_mass <- pars$var[, 'slurry_mass']

  # Get timing of intervals
  timelist <- makeTimeList(pars, times, days, delta_t)
  n_int <- length(timelist)
  
  # Empty data frame for holding results
  dat <- NULL

  # Time trackers
  # Time remaining to run
  t_rem <- days
  # Time that has already run
  t_run <- 0

  # NTS: below comments are no longer completely true
  # Note: Removals, slurry_mass, wash_water, slurry_prod_rate_t, and t_ints all have same length,
  # Note: but all except _mass have placeholder in first position
  
  # Start the time (emptying) loop
  for (i in 2:n_int) {

    # Sort out call duration
    t_call <- min(max(timelist[[i]]), t_rem)

    # Fill in current pars from var
    # Also adds 2 temperature-dependent derivative vectors
    pars <- updateVarPars(pars, y, i - 1)

    # Calculate log Ka for speciation
    pars <- calcKa(pars)
    
    # Create empty (0) y.eff vector because washing could occur, and dat needs columns
    y.eff <- 0 * emptyStore(y)$eff

    # If there is a removal event, remove slurry before calling up ODE solver
    if (pars$removal) {
      y <- emptyStore(y, resid_mass = pars$resid_mass, resid_enrich = pars$resid_enrich)
      y.eff <- y$eff
      y <- y$store
    }

    # Get times for lsoda() call
    # Need some care with times to make sure t_call is last one in case it is not multiple of delta_t
    tt <- timelist[[i]]

    # Add run time to pars so rates() can use actual time to calculate temp_C and pH
    pars$t_run <- t_run
    
    # Call up ODE solver
    out <- deSolve::lsoda(y = y, 
                       times = tt, 
                       rates, 
                       parms = pars)
     
    # Change format of output and drop first (time 0) row (duplicated in last row of previous)
    if (i == 2) {
      out <- data.frame(out)
    } else {
      out <- data.frame(out[-1, , drop = FALSE])
    }

    # Extract new state variable vector from last row of lsoda output
    y <- getLastState(out, y)

    # Add effluent results
    out[, names(y.eff)] <- 0
    out[nrow(out), names(y.eff)] <- y.eff
 
    # Clean up and stack output with earlier results
    dat <- addOut(dat, out)
    
    # Update time remaining and total time run so far
    t_rem <- t_rem - t_call
    t_run <- t_run + t_call
    
  }

  # Drop times from slurry mass data that were not requested in output
  if (!is.null(times)) {
    dat <- dat[dat$time %in% c(times, days), ]
  }

  return(dat)
}

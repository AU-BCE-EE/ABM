abm_variable <-
  function(days, 
           delta_t, 
           times, 
           y, 
           pars, 
           warn, 
           temp_C_fun = temp_C_fun, 
           pH_fun = pH_fun, 
           slurry_mass_approx) {
    
  pars$abm_regular <- FALSE
  
  # Some warnings about unused inputs
  if (!is.na(pars$wash_water) & pars$wash_water != 0) {
    warning('Fixed wash_water value of ', pars$wash_water, '\nwill be ignored because variable slurry input is used.')
  }

  pars <- fixSlurryMass(pars, days, slurry_mass_approx)
  removals <- pars$removals

  # Extract washing mass
  wash_water <- getWashWater(pars)
  
  # Extract slurry_mass vector for use in emptying calculations
  slurry_mass <- pars$slurry_mass[, 'slurry_mass']

  # Determine variable slurry production rate for each time interval
  slurry_prod_rate_t <- getSlurryProd(pars)

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

  # Note: Removals, slurry_mass, wash_water, slurry_prod_rate_t, and t_ints all have same length,
  # Note: but all except _mass have placeholder in first position
  # Start the time (emptying) loop

  for (i in 2:n_int) {

    # Sort out call duration
    t_call <- min(max(timelist[[i]]), t_rem)

    # Fill in slurry_prod_rate
    pars$slurry_prod_rate <- slurry_prod_rate_t[i]
    
    # Create empty (0) y.eff vector because washing could occur, and dat needs columns
    yy <- emptyStore(y, resid_mass = 0, resid_enrich = 0)
    y.eff <- 0 * yy$eff

    # If there is a removal event, remove slurry before calling up ODE solver
    if (removals[i]) {
      if (slurry_mass_approx == 'late') {
        j <- i - 1
      } else {
        j <- i
      }
      y <- emptyStore(y, resid_mass = slurry_mass[j], resid_enrich = pars$resid_enrich)
      y.eff <- y$eff
      y <- y$store

      # Washing, increase slurry mass
      # Occurs after emptying here 
      if (wash_water[i] > 0) {
        y['slurry_mass'] <- y['slurry_mass'] + wash_water[i]

        # And empty again, leaving same residual mass as always, and enriching for particles
        y <- emptyStore(y, resid_mass = slurry_mass[j], resid_enrich = pars$resid_enrich)
        y.eff <- y$eff
        y <- y$store
      }
    }

    # This might not be possible with above code
    if (t_call <= 0) stop('t_call < 0. slk409.')

    # Get times for lsoda() call
    # Need some care with times to make sure t_call is last one in case it is not multiple of delta_t
    tt <- timelist[[i]]

    # Add run time to pars so rates() can use actual time to calculate temp_C and pH
    pars$t_run <- t_run
    
    # Call up ODE solver
    out <- deSolve::lsoda(y = y, 
                       times = tt, 
                       rates, 
                       parms = pars, 
                       temp_C_fun = temp_C_fun, 
                       pH_fun = pH_fun)
     
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
    dat <- addOut(dat, out, t_add = t_run, y.eff = y.eff)
    
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

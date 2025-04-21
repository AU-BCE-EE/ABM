abm_variable <-
  function(days, delta_t, times, y, pars, warn, temp_C_fun = temp_C_fun, pH_fun = pH_fun, 
           SO4_inhibition_fun = SO4_inhibition_fun, conc_fresh_fun = conc_fresh_fun, xa_fresh_fun = xa_fresh_fun, slurry_mass_approx) {
    
  #initialize dat for storage of results to speed up bind_rows
  dat <- as.data.frame(matrix(NA, nrow = days * 2, ncol = 400)) # slow speed, but much faster than before
  
  pars$abm_regular <- FALSE
  
  # Some warnings about unused inputs
  if (!is.na(pars$wash_water) & pars$wash_water != 0) {
    warning('Fixed wash_water value of ', pars$wash_water, '\nwill be ignored because variable slurry input is used.')
  }

  # Cannot have no slurry present because is used in all concentration calculations (NTS: could change this)
  pars$slurry_mass[pars$slurry_mass[, 'slurry_mass'] == 0, 'slurry_mass'] <- 1E-10

  # Trim unused times
  pars$slurry_mass <- pars$slurry_mass[pars$slurry_mass$time <= days, ]

  # Check for sorted time
  if (is.unsorted(pars$slurry_mass$time)) {
    stop('Column `time` must be sorted when `slurry_mass` is time-dependent, but it is not: ',
         pars$slurry_mass$time)
  }
  
  # If simulation continues past slurry_mass data frame time, set following slurry production (addition) to zero
  # Extend last slurry mass all the way 
  if (pars$slurry_mass[nrow(pars$slurry_mass), 'time'] < days) {
    t_end <- days
    pars$slurry_mass <- rbind(pars$slurry_mass, pars$slurry_mass[nrow(pars$slurry_mass), ])
    pars$slurry_mass[nrow(pars$slurry_mass), 'time'] <- days
    # But make sure washing is not repeated!
    if (ncol(pars$slurry_mass) > 2) {
      pars$slurry_mass[nrow(pars$slurry_mass), 3:ncol(pars$slurry_mass)] <- 0
    }
  }

  # Sort out times returned by ODE solver
  if (is.null(times)) {
    times <- seq(0, days, by = delta_t)
  }

  # Include days argument in times vector
  times <- sort(unique(c(times, days)))
  times.orig <- times
  times <- sort(unique(c(times, pars$slurry_mass$time)))
  droptimes <- times[!times %in% times.orig]

  # Adjust slurry_mass for rain and evaporation, so that it is actually the mass of *slurry* added or removed (excluding downstream addition or removal of water) 
  # First get net precip/evap addition as a daily rate in kg/d (mm/d * m2 * density = kg/d)
  pe_netr <- (pars$rain - pars$evap) * pars$area * (pars$dens / 1000) 
  slurry_mass_orig <- pars$slurry_mass$slurry_mass
  pars$slurry_mass$slurry_mass <- pars$slurry_mass$slurry_mass - cumsum(c(0, diff(pars$slurry_mass$time) * pe_netr))
  if (any(pars$slurry_mass$slurry_mass < 0)) {
    stop('Negative slurry mass values after adjustment for rain and evaporation.\n  Check parameters and try again.')
  }
  
  # Determine slurry removal quantity in each time interval
  # Note final 0--alignment is a bit tricky
  if (slurry_mass_approx == 'late') {
    removals <- - c(0, 0, diff(pars$slurry_mass[-nrow(pars$slurry_mass), 'slurry_mass'])) > 0
  } else if (slurry_mass_approx == 'early') {
    removals <- - c(0, diff(pars$slurry_mass[, 'slurry_mass'])) > 0
  } else if (slurry_mass_approx == 'mid') {
    # Get midpoint time
    ir <- which(- c(0, diff(pars$slurry_mass[, 'slurry_mass'])) > 0)
    tt <- (pars$slurry_mass[ir, 'time']  + pars$slurry_mass[ir - 1, 'time']) / 2
    # And copy slurry mass
    mm <- pars$slurry_mass[ir, 'slurry_mass']
    pars$slurry_mass <- rbind(pars$slurry_mass, data.frame(time = tt, slurry_mass = mm))
    pars$slurry_mass <- pars$slurry_mass[order(pars$slurry_mass$time), ]
    # Then apply early (late??) removal method
    removals <- - c(0, 0, diff(pars$slurry_mass[-nrow(pars$slurry_mass), 'slurry_mass'])) > 0
    #removals <- - c(0, diff(pars$slurry_mass[, 'slurry_mass'])) > 0
    slurry_mass_approx <- 'late'
  } 

  # Extract slurry_mass for use in emptying calculations
  slurry_mass <- pars$slurry_mass[, 'slurry_mass']
  
  # For removal/emptying events, which are instant, set slurry_mass back to original
  if (slurry_mass_approx %in% c('late', 'mid')) {
    slurry_mass[c(removals[-1], FALSE)] <- slurry_mass_orig[c(removals[-1], FALSE)]
  } else if (slurry_mass_approx == 'early') {
    slurry_mass[removals] <- slurry_mass_orig[removals]
  } 
  
  # Extract washing mass
  if ('wash_water' %in% names(pars$slurry_mass)) {
    wash_water <- pars$slurry_mass[, 'wash_water']
  } else {
    wash_water <- 0 * slurry_mass
  }

  # NTS: can put evap here, after calculating daily temperature
  # NTS: call evap function, check old commit, see wind speed branch
  # wash_water <- NULL
  # ifelse(length(pars$wash_water) <= 1, wash_water <- 0, wash_water <- pars$wash_water[, 'wash_water'])

  # Determine variable slurry production rate for each time interval
  slurry_prod_rate_t <- c(0, diff(pars$slurry_mass[, 'slurry_mass']) / diff(pars$slurry_mass[, 'time'])) 
  slurry_prod_rate_t[slurry_prod_rate_t < 0] <- 0
  slurry_prod_rate_t[!is.finite(slurry_prod_rate_t)] <- 0

  # Get number of timing of intervals
  # Note about time: 1) All simulations start at 0, 2) days must be at least as long as mass data
  # Note that this works even with t_end = NULL (is ignored)
  # Note the "dummy" placeholder in position 1 (and extra + 1 in n_int)
  n_int <- nrow(pars$slurry_mass)
  st <- pars$slurry_mass[, 'time']
  timelist <- cumtime <- as.list(0 * removals)
  for (i in 2:n_int) {
    tt <- times[times > st[i - 1] & times <= st[i]] # slow speed of this line
    # Simulation intervals should end exactly at emptying time, so it is added here through st[i] for cases where there is not alignment
    tt <- unique(c(tt, st[i]))
    if (length(tt) == 0) { 
      # Not clear when this might happen, but it isn't good
      stop('No times for interval (row) ', i, ' in time list for some reason. hatty917')
    } else {
      cumtime[[i]] <- max(tt)
      tt <- tt - cumtime[[i - 1]]
      tt <- unique(c(0, tt))
    }
    timelist[[i]] <- tt
  }
  
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

    # Get number of microbial groups
    n_mic <- length(pars$grps)

    # Fill in slurry_prod_rate
    pars$slurry_prod_rate <- slurry_prod_rate_t[i]

    # With no removal
    # Call lsoda and advance if time changes (i.e., if there is not a removal)
    # Note that t_call should always by 0 when there is a removal, so skipping this code is correct

    # If there is a removal event, remove slurry before calling up ODE solver
    if (removals[i]) {
      if (slurry_mass_approx == 'late') {
        j <- i - 1
      } else {
        j <- i
      }
      y <- emptyStore(y, resid_mass = slurry_mass[j], resid_enrich = pars$resid_enrich)
      y.eff <- y[grepl('_eff$', names(y))]
      y <- y[!grepl('_eff$', names(y))]

      # Washing, increase slurry mass
      # Typically occurs after emptying here (wash_int ignored)
      # But could take place at the *end* of *any* interval
      if (wash_water[i] > 0) {
        y['slurry_mass'] <- y['slurry_mass'] + wash_water[i]

        # And empty again, leaving same residual mass as always, and enriching for particles
        y <- emptyStore(y, resid_mass = slurry_mass[j], resid_enrich = pars$resid_enrich)
        y.eff <- y.eff + y[grepl('_eff$', names(y))]
        y <- y[!grepl('_eff$', names(y))]
        out[nrow(out), names(y)] <- y
      }
    }

    # This might not be possible with above code
    if (t_call <= 0) stop('t_call < 0. slkj4098sh, check if "time" is duplicated in slurry mass data frame')

    # Need some care with times to make sure t_call is last one in case it is not multiple of delta_t
    times <- timelist[[i]]

    # Add run time to pars so rates() can use actual time to calculate temp_C and pH
    pars$t_run <- t_run

    # Call up ODE solver
    out <- deSolve::lsoda(y = y, times = times, rates, parms = pars, temp_C_fun = temp_C_fun, 
                          pH_fun = pH_fun, SO4_inhibition_fun = SO4_inhibition_fun, 
                          conc_fresh_fun = conc_fresh_fun, xa_fresh_fun = xa_fresh_fun)
    
    # Change format of output and drop first (time 0) row (duplicated in last row of previous)
    if (i == 2) {
      out <- data.frame(out)
    } else {
      out <- data.frame(out[-1, , drop = FALSE])
    }

    # Fix some names (NTS: can we sort this out in rates()? I have tried and failed)
    names(out)[1:length(pars$qhat_opt) + 1] <- names(pars$qhat_opt)

    # Extract new state variable vector from last row of lsoda output
    y <- unlist(out[nrow(out), 1:(length(c(pars$qhat_opt)) + 32) + 1])

    # Change time in output to cumulative time for complete simulation
    out$time <- out$time + t_run
    
    # Create empty (0) y.eff vector because washing could occur, and dat needs columns
    yy <- emptyStore(y, resid_mass = 0, resid_enrich = 0)
    y.eff <- 0 * yy[grepl('_eff$', names(yy))]

    # Add effluent results
    out[, names(y.eff)] <- 0
    out[nrow(out), names(y.eff)] <- y.eff
 
    # Add results to earlier ones
    dat <- dplyr::bind_rows(dat, out)
    #dat <- rbind(dat, out)
    
    # Update time remaining and total time run so far
    t_rem <- t_rem - t_call
    t_run <- t_run + t_call

  }

  # Drop times from slurry mass data that were not requested in output

  dat <- dat[!dat$time %in% droptimes, ]

  return(dat)
}

# Simulation with slurry mass given as a data frame (time-variable filling and emptying)
# Not exported, called by abm()

abm_variable <- function(days, delta_t, times, y, pars, warn, slurry_mass_approx) {
    
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
  extended <- pars$slurry_mass[nrow(pars$slurry_mass), 'time'] < days
  if (extended) {
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
  # Cannot combine this with mid approach--technically impossible and practically difficult
  if (slurry_mass_approx != 'mid') {
    pe_netr <- (pars$rain - pars$evap) * pars$area * (pars$dens / 1000)
    slurry_mass_orig <- pars$slurry_mass$slurry_mass
    pars$slurry_mass$slurry_mass <- pars$slurry_mass$slurry_mass - cumsum(c(0, diff(pars$slurry_mass$time) * pe_netr))
  } else {
    if (any(c(pars$rain != 0, pars$evap != 0))) {
      warning('Setting rain and evap to 0 because "mid" method was selected.')
    }
    pars$rain <- pars$evap <- 0
    pe_netr <- 0
  }

  if (any(pars$slurry_mass$slurry_mass < 0)) {
    warning('Negative slurry mass values after adjustment for rain and evaporation.\n  Check parameters and try again.')
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
    
    # if wash_water exists for some reason, we need to add it when rbinding.
    if ("wash_water" %in% colnames(pars$slurry_mass)) {
      if (length(mm) == 0) {
        wash_water_fill <- mm
      } else {
        wash_water_fill <- 0
      }
      pars$slurry_mass <- rbind(pars$slurry_mass, data.frame(time = tt, slurry_mass = mm, wash_water = wash_water_fill))
    } else {
      pars$slurry_mass <- rbind(pars$slurry_mass, data.frame(time = tt, slurry_mass = mm))
    }
    
    pars$slurry_mass <- pars$slurry_mass[order(pars$slurry_mass$time), ]
    # Then apply late removal method
    removals <- - c(0, 0, diff(pars$slurry_mass[-nrow(pars$slurry_mass), 'slurry_mass'])) > 0
  } else {
    stop('approx_method["slurry_mass"] must be "early", "late", or "mid", but is "', slurry_mass_approx, '".')
  }

  # Extract slurry_mass for use in emptying calculations
  slurry_mass <- pars$slurry_mass[, 'slurry_mass']

  # For removal/emptying events, which are instant, set slurry_mass back to original
  if (slurry_mass_approx == 'late') {
    slurry_mass[c(removals[-1], FALSE)] <- slurry_mass_orig[c(removals[-1], FALSE)]
  } else if (slurry_mass_approx == 'early') {
    slurry_mass[removals] <- slurry_mass_orig[removals]
  } else if (slurry_mass_approx == 'mid') {
    slurry_mass_approx <- 'late'
  }

  # Extract washing mass
  if ('wash_water' %in% names(pars$slurry_mass)) {
    wash_water <- pars$slurry_mass[, 'wash_water']
  } else {
    wash_water <- 0 * slurry_mass
  }

  # Washing only takes place at a removal (emptying) event (see loop below), so wash_water anywhere else would be
  # ignored without notice. Warn about it.
  if (warn) {
    lost <- which(!is.na(wash_water) & wash_water > 0 & !removals)
    if (length(lost) > 0) {
      lt <- signif(pars$slurry_mass[lost, 'time'], 6)
      where <- if (slurry_mass_approx == 'early') 'the row where slurry_mass drops' else 'the row after the drop in slurry_mass'
      warning(length(lost), ' wash_water value(s) in slurry_mass not applied (', signif(sum(wash_water[lost]), 4), ' kg in total, at time ',
              paste(head(lt, 10), collapse = ', '), if (length(lt) > 10) ', ...', ').\n',
              '  Washing only takes place at a slurry removal event: wash_water must sit on a removal, i.e. on ', where,
              ' (approx_method slurry_mass = "', slurry_mass_approx, '").\n',
              '  A drop only counts as a removal if it is larger than net evaporation minus rain since the previous row, because\n',
              '  slurry_mass is corrected for these. At the wash row, give slurry_mass as the level after emptying (without the\n',
              '  wash water); otherwise a rise from wash water is read as fresh slurry.')
    }
  }

  # Determine variable slurry production rate for each time interval
  slurry_prod_rate_t <- c(0, diff(pars$slurry_mass[, 'slurry_mass']) / diff(pars$slurry_mass[, 'time']))

  # Rain for each time interval
  # In filling intervals where the given slurry mass rises less than net rain - evaporation, slurry addition would be negative
  # Then set slurry addition to 0 and reduce rain so simulated slurry mass follows the given slurry mass
  rain_t <- rep(pars$rain, length(slurry_prod_rate_t))
  dt_sm <- c(NA, diff(pars$slurry_mass[, 'time']))
  # Given (unadjusted) mass increase rate
  dm_given <- c(NA, diff(pars$slurry_mass[, 'slurry_mass'])) / dt_sm + pe_netr
  short <- which(is.finite(slurry_prod_rate_t) & slurry_prod_rate_t < 0 & is.finite(dm_given) & dm_given >= 0)
  # Not for the added final interval when simulation continues past given slurry mass data (no data there)
  if (extended) short <- short[short != length(slurry_prod_rate_t)]
  if (length(short) > 0) {
    rain_t[short] <- dm_given[short] / (pars$area * pars$dens / 1000) + pars$evap
    if (warn) {
      st_sm <- pars$slurry_mass[, 'time']
      excess <- sum((pars$rain - rain_t[short]) * pars$area * pars$dens / 1000 * dt_sm[short])
      warning('In ', length(short), ' interval(s) of slurry_mass (e.g., time ', st_sm[short[1] - 1], ' to ', st_sm[short[1]],
              ') the given slurry mass increases less than net rain - evaporation (', signif(pe_netr, 3), ' kg/d).\n',
              '  Slurry addition was set to 0 and rain reduced in these intervals so simulated slurry mass follows the given values\n',
              '  (', signif(excess, 3), ' kg less rain in total). Check rain and evap.')
    }
  }

  slurry_prod_rate_t[slurry_prod_rate_t < 0] <- 0
  slurry_prod_rate_t[!is.finite(slurry_prod_rate_t)] <- 0

  # Get number of timing of intervals
  # Note about time: 1) All simulations start at 0, 2) days must be at least as long as mass data
  # Note the "dummy" placeholder in position 1 (and extra + 1 in n_int)
  n_int <- nrow(pars$slurry_mass)
  st <- pars$slurry_mass[, 'time']
  # Output times within each interval (st[i - 1], st[i]], found in one step
  # Simulation intervals should end exactly at emptying time, so st[i] is added for cases where there is not alignment
  tgrp <- split(times, factor(findInterval(times, st, left.open = TRUE) + 1, levels = 1:n_int))
  timelist <- vector('list', n_int)
  cumtime <- c(0, st[-1])
  for (i in 2:n_int) {
    tt <- unique(c(tgrp[[i]], st[i]))
    tt <- tt - cumtime[i - 1]
    timelist[[i]] <- unique(c(0, tt))
  }
  
  # Number of state variables
  n_state <- length(y)
  n_mic <- pars$n_mic

  # Parameter indices, output names etc. for rates_cpp()
  pars <- prepRatesPars(pars)

  # Names of effluent columns (all state variables except cumulative ones)
  eff_names <- names(emptyStore(y, resid_mass = 0, resid_enrich = c(OM = 0, mic = 0)))
  eff_names <- eff_names[grepl('_eff$', eff_names)]

  # List for holding results (combined at end, much faster than growing a data frame)
  outs <- vector('list', n_int)
  last_out <- NULL
  pending_eff <- NULL

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

    # Fill in slurry_prod_rate and rain
    pars$slurry_prod_rate <- slurry_prod_rate_t[i]
    pars$rain <- rain_t[i]

    # If there is a removal event, remove slurry before calling up ODE solver
    if (removals[i]) {
      
      if (slurry_mass_approx == 'late') {
        j <- i - 1
      } else {
        j <- i
      }
      
      y <- emptyStore(y, resid_mass = slurry_mass[j], resid_enrich = pars$resid_enrich)
      which_eff <- grepl('_eff$', names(y))
      y.eff <- y[which_eff]
      y <- y[!which_eff]

      # Washing, increase slurry mass
      # Typically occurs after emptying here (wash_int ignored)
      # But could take place at the *end* of *any* interval
      if (wash_water[i] > 0) {
        y['slurry_mass'] <- y['slurry_mass'] + wash_water[i]
        
        # And empty again, leaving same residual mass as always, and enriching for particles
        y <- emptyStore(y, resid_mass = slurry_mass[j], resid_enrich = pars$resid_enrich)
        which_eff <- grepl('_eff$', names(y))
        y.eff <- y.eff + y[which_eff]
        y <- y[!which_eff]
      }

      # Track removed COD, C, and N in cumulative state variables (for mass balance)
      y <- addEffCum(y, y.eff, pars)

      # Record effluent in last output row before removal (or first row if there is no earlier output)
      if (is.null(last_out)) {
        pending_eff <- y.eff
      } else {
        outs[[last_out]][nrow(outs[[last_out]]), names(y.eff)] <- y.eff
      }
    }

    # This might not be possible with above code
    if (t_call <= 0) stop('t_call < 0. slkj4098sh, check if "time" is duplicated in slurry mass data frame')

    times <- timelist[[i]]
    
    # Add run time to pars so rates() can use actual time to calculate temp_C and pH
    pars$t_run <- t_run

    # Call up ODE solver
    out <- runLsoda(y = y, times = times, pars = pars)
    out <- matrix(out, nrow = nrow(out), dimnames = dimnames(out))

    # Extract new state variable vector from last row of lsoda output
    y <- out[nrow(out), 1:n_state + 1]

    # Microbial group names (instead of e.g., xa.m1)
    colnames(out)[1:n_mic + 1] <- names(pars$qhat_opt)

    # Drop first (time 0) row (duplicated in last row of previous)
    if (!is.null(last_out)) {
      out <- out[-1, , drop = FALSE]
    }

    # Change time in output to cumulative time for complete simulation
    out[, 'time'] <- out[, 'time'] + t_run
    
    # Add effluent columns (filled in at removal events)
    eff <- matrix(0, nrow = nrow(out), ncol = length(eff_names), dimnames = list(NULL, eff_names))
    if (!is.null(pending_eff)) {
      eff[1, names(pending_eff)] <- pending_eff
      pending_eff <- NULL
    }
    outs[[i]] <- cbind(out, eff)
    last_out <- i
    
    # Update time remaining and total time run so far
    t_rem <- t_rem - t_call
    t_run <- t_run + t_call

  }

  dat <- as.data.frame(do.call(rbind, outs))

  # Drop times from slurry mass data that were not requested in output
  # Effluent from emptying at times that are not returned is added to the next returned time
  dat <- filterTimes(dat, keep = !dat$time %in% droptimes)

  return(dat)
}

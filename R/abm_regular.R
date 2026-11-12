# Simulation with fixed slurry production rate and regular emptying (and optional washing)
# Not exported, called by abm()

abm_regular <- function(days, delta_t, times_regular, y, pars, starting = NULL) { 
  
  empty_int <- pars$empty_int
  
  # if empty interval is set to 0 or NA it means the storage is never emptied. 
  if (is.na(empty_int) || empty_int == 0) empty_int <- days + 1
  
  # Figure out time intervals for loop
  if (!is.na(pars$wash_int) && pars$wash_water > 0) {  
    wash_int <- pars$wash_int
    rest_d <- pars$rest_d
  } else {
    wash_int <- Inf
    rest_d <- 0
  }

  # Continue sorting out intervals
  i <- 0
  t_int <- 0
  t_nowash <- 0
  wash <- FALSE

  # Each interval is either 1) the fixed empty_int or if time between washings would be exceeded, 2) time to get to a washing event
  while (sum(t_int, wash * rest_d) < days) {
    i <- i + 1
    t_int[i] <- min(wash_int - t_nowash, empty_int, days - sum(t_int, wash * rest_d))
    if (t_int[i] == wash_int - t_nowash) {
      wash[i] <- TRUE
      t_nowash <- 0
    } else {
      wash[i] <- FALSE
      t_nowash <- t_nowash + t_int[i]
    }
  }

  # Number of empty or wash intervals
  n_int <- length(t_int)

  # Number of state variables
  n_state <- length(y)
  n_mic <- pars$n_mic

  # Parameter indices, output names etc. for rates_cpp()
  pars <- prepRatesPars(pars)

  # List for holding results (combined at end, much faster than growing a data frame)
  outs <- vector('list', n_int)

  # Time trackers
  # Time remaining to run
  t_rem <- days
  # Time that has already run
  t_run <- 0
  
  # Start the time (emptying) loop
  for (i in 1:n_int) {
    
    # Sort out call duration
    t_call <- t_int[i]
    
    # Need some care with times to make sure t_call is last one in case it is not multiple of delta_t
    # Also include any requested output times (times_regular) within this interval
    tr <- times_regular[times_regular > t_run & times_regular < t_run + t_call] - t_run
    times <- sort(unique(round(c(seq(0, t_call, by = min(t_rem, delta_t)), t_call, tr), 5)))
    
    # Add run time to pars so rates() can use actual time to calculate temp_C and pH
    pars$t_run <- t_run

    # Call up ODE solver
    out <- deSolve::lsoda(y = y, times = times, rates_cpp, parms = pars, rtol = 1E-5, atol = 1E-5)
    out <- matrix(out, nrow = nrow(out), dimnames = dimnames(out))
    
    # Extract new state variable vector from last row
    y <- out[nrow(out), 1:n_state + 1]
    
    # Empty channel (instantaneous changes at end of day) in preparation for next lsoda call
    y <- emptyStore(y, resid_mass = pars$resid_mass, resid_enrich = pars$resid_enrich)
    which_eff <- grepl('_eff$', names(y))
    y.eff <- y[which_eff]
    y <- y[!which_eff]
    
    # Washing, increase slurry mass
    if (wash[i]) {
      y['slurry_mass'] <- y['slurry_mass'] + pars$wash_water

      # And empty again, leaving same residual mass as always, and enriching for particles (same as abm_variable())
      y <- emptyStore(y, resid_mass = pars$resid_mass, resid_enrich = pars$resid_enrich)
      which_eff <- grepl('_eff$', names(y))
      y.eff <- y.eff + y[which_eff]
      y <- y[!which_eff]
    }

    # Add effluent results (in last row, i.e., at time of emptying, before removal)
    eff <- matrix(0, nrow = nrow(out), ncol = length(y.eff), dimnames = list(NULL, names(y.eff)))
    eff[nrow(out), ] <- y.eff
    out <- cbind(out, eff)

    # Run for rest period with no influent after washing (but not past end of simulation)
    rest_i <- if (wash[i]) min(rest_d, days - t_run - t_call) else 0
    if (rest_i > 0) {
      t_run_r <- t_run + t_call
      tr <- times_regular[times_regular > t_run_r & times_regular < t_run_r + rest_i] - t_run_r
      times <- sort(unique(round(c(seq(0, rest_i, by = delta_t), rest_i, tr), 5)))
      parsr <- pars
      parsr$slurry_prod_rate <- 0
      # Actual time, for temperature etc.
      parsr$t_run <- t_run_r
      outr <- deSolve::lsoda(y = y, times = times, rates_cpp, parms = parsr, rtol = 1E-5, atol = 1E-5)
      outr <- matrix(outr, nrow = nrow(outr), dimnames = dimnames(outr))
      # Extract new state variable vector from last row
      y <- outr[nrow(outr), 1:n_state + 1]
      # Correct time in outr and combine with main output (no effluent)
      outr[, 'time'] <- outr[, 'time'] + t_call
      effr <- matrix(0, nrow = nrow(outr), ncol = ncol(eff), dimnames = list(NULL, colnames(eff)))
      out <- rbind(out, cbind(outr, effr))
    }
    
    # Microbial group names (instead of e.g., xa.m1)
    colnames(out)[1:n_mic + 1] <- names(pars$qhat_opt)

    # Change time in output to cumulative time for complete simulation
    out[, 'time'] <- out[, 'time'] + t_run
  
    outs[[i]] <- out
    
    # Update time remaining and time run so far
    t_rem <- t_rem - t_call - rest_i
    t_run <- t_run + t_call + rest_i

  }

  dat <- as.data.frame(do.call(rbind, outs))

  # Return only requested times, if given
  # Effluent from emptying at times that are not returned is added to the next returned time
  if (!is.null(times_regular)) {
    times_regular <- sort(unique(c(times_regular, days)))
    dat <- filterTimes(dat, keep = round(dat$time, 5) %in% round(times_regular, 5))
  }

  return(dat)

}

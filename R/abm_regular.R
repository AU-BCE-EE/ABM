abm_regular <- function(days, delta_t, times_regular, y, pars, starting = NULL, temp_C_fun = temp_C_fun, pH_fun = pH_fun, 
                        conc_fresh_fun = conc_fresh_fun, xa_fresh_fun = xa_fresh_fun) { 
  
  #initialize dat for storage of results to speed up bind_rows
  dat <- as.data.frame(matrix(NA, nrow = days * 2, ncol = 400)) 
  
  empty_int <- pars$empty_int
  
  # if empty interval is set to 0 or NA it means the storage is never emptied. 
  if(empty_int == 0 || is.na(empty_int)) empty_int <- days + 1
  
  if (pars$slurry_mass == 0) pars$slurry_mass <- 1E-10

  # Figure out time intervals for loop
  if (!is.na(pars$wash_int) && pars$wash_water > 0) {  
    wash_int <- pars$wash_int
    rest_d <- pars$rest_d
  } else {
    wash_int <- Inf
    rest_d <- 0
  }
  wash_rest_int <- wash_int + rest_d

  # Continue sorting out intervals
  i <- 0
  t_int <- 0
  t_nowash <- 0
  wash <- FALSE

  # Continute . . .
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

  # Empty data frame for holding results
  dat <- NULL

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
    times <- sort(unique(round(c(seq(0, t_call, by = min(t_rem, delta_t)), t_call), 5)))

    # Add run time to pars so rates() can use actual time to calculate temp_C and pH
    pars$t_run <- t_run
    pars$t_call <- t_call
    pars$times <- times
    
    pars$p_idx <- pars_indices(pars)
    pars$temp_C_fun <- temp_C_fun
    pars$pH_fun <- pH_fun
    pars$CTM_cpp <- CTM_cpp
    pars$H2SO4_titrate <- H2SO4_titrat
    pars$xa_fresh_fun <- xa_fresh_fun
    pars$conc_fresh_fun <- conc_fresh_fun
    # Call up ODE solver
    #cat(t_rem, '\n')

    out <- deSolve::lsoda(y = y, times = times, rates_cpp, parms = pars,rtol = 1E-5, atol = 1E-5)

    # Get number of microbial groups
    n_mic <- length(pars$n_mic)

    # Extract new state variable vector from last row
    y <- out[nrow(out), 1:(length(c(pars$qhat_opt)) + 27) + 1]

    
    # Empty channel (instantaneous changes at end of day) in preparation for next lsoda call
    y <- emptyStore(y, resid_mass = pars$resid_mass, resid_enrich = pars$resid_enrich)
    y.eff <- y[grepl('_eff$', names(y))]
    y <- y[!grepl('_eff$', names(y))]
    #y[which.effluent] <- (y.before - y)[gsub('_eff$', '', names(y)[which.effluent])]
   
    # Washing, increase slurry mass
    if (wash[i]) {
      y['slurry_mass'] <- y['slurry_mass'] + pars$wash_water

      # And empty again, leaving same residual mass as always, and enriching for particles
      y <- emptyStore(y, resid_mass = pars$resid_mass, resid_enrich = pars$resid_enrich)
      y.eff <- y.eff + y[grepl('_eff$', names(y))]
      y <- y[!grepl('_eff$', names(y))]

      # Run for rest period with no influent 
      if (pars$rest_d > 0) {
        times <- seq(0, pars$rest_d, delta_t)
        parsr <- pars
        parsr$slurry_prod_rate <- 0
        outr <- deSolve::lsoda(y = y, times = times, rates_cpp, parms = parsr, rtol = 1E-5, atol = 1E-5)
        # Extract new state variable vector from last row
        y <- outr[nrow(outr), 1:(length(c(pars$qhat_opt)) + 27) + 1]
        # Correct time in outr and combine with main output
        outr[, 'time'] <- outr[, 'time'] + out[nrow(out), 'time']
        out <- rbind(out, outr)
      }
    }
    
    # Change format of output
    # Do not drop first (time 0) row
    out <- data.frame(out)

    # Add effluent results
    out[, names(y.eff)] <- 0
    out[nrow(out), names(y.eff)] <- y.eff

    # Fix some names (NTS: can we sort this out in rates()? I have tried and failed)
    names(out)[1:length(pars$qhat_opt) + 1] <- names(pars$qhat_opt)

    # Change time in output to cumulative time for complete simulation
    out$time <- out$time + t_run
  
    # Add results to earlier ones
    dat <- dplyr::bind_rows(dat, out)
    
    # Update time remaining and time run so far
    t_rem <- t_rem - t_call - wash[i] * rest_d
    t_run <- t_run + t_call + wash[i] * rest_d

  }

  if (!is.null(times_regular)) {
    times_regular <- sort(unique(c(times_regular, days)))
    dat <- dat[dat$time %in% times_regular, ]
  }

  return(dat)

}

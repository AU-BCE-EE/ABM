abm_regular <- function(days, delta_t, y, pars, conc_fresh_sel_fun, temp_fun, pH_fun) {

  pars$abm_regular <- TRUE
  
  empty_int <- pars$empty_int

  # Cannot have no slurry present because is used in all concentration calculations (NTS: could change this)
  if (pars$slurry_mass == 0) {
    pars$slurry_mass <- 1E-10
  }

  # Sorting out intervals
  t_int <- c(rep(empty_int, days %/% empty_int), days %% empty_int)

  # Number of empty intervals
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

    # Add run time to pars so rates() can use actual time to calculate temp and pH
    pars$t_run <- t_run

    # Call up ODE solver
    #cat(t_rem, '\n')
    out <- deSolve::lsoda(y = y, times = times, rates, parms = pars, conc_fresh_sel_fun = conc_fresh_sel_fun, temp_fun = temp_fun, pH_fun = pH_fun)

    # Get number of microbial groups
    n_mic <- length(pars$qhat_opt)

    # Extract new state variable vector from last row
    y <- out[nrow(out), 1:(length(y)) + 1]

    # Empty channel (instantaneous changes at end of day) in preparation for next lsoda call
    if (t_call >= empty_int) {
      if (y['slurry_mass'] > pars$resid_mass) {
        resid_frac <- pars$resid_mass / y['slurry_mass']
        resid_xa <- logistic(logit(resid_frac) + pars$resid_enrich)

        y['slurry_mass'] <- max(resid_frac * y['slurry_mass'], 1E-10)

        if (pars$mix) {
          # Note invert argument below
          empty_names <- grep('^xa\\.|cum|slurry_mass', names(y), value = TRUE, invert = TRUE)
          # If lagoon contents are mixed microbial biomass can still be retained beyond resid_frac
          # If this isn't desired, set resid_enrich = 0
          resid_xa <- logistic(logit(resid_frac) + pars$resid_enrich)
          y[1:n_mic] <- resid_xa * y[1:n_mic]
        } else {
          # Note invert argument below
          empty_names <- grep('^xa\\.|cum|slurry_mass|sed$', names(y), value = TRUE, invert = TRUE)
          # Without mixing, microbial bioass is assumed to be completely sedimented and is completely retained
          # All sediment is too of course
        }
        y[empty_names] <- sapply(y[empty_names], function(x) resid_frac * x) # apply resid_frac to other state variables

      } else {
        warning('Emptying skipped because of low slurry level.')
      }
    }
   
    # Change format of output
    # Do not drop first (time 0) row
    out <- data.frame(out)

    # Fix some names (NTS: can we sort this out in rates()? I have tried and failed)
    names(out)[1:length(pars$qhat_opt) + 1] <- names(pars$qhat_opt)

    # Change time in output to cumulative time for complete simulation
    out$time <- out$time + t_run
  
    # Add results to earlier ones
    dat <- rbind(dat, out)
 
    # Update time remaining and time run so far
    t_rem <- t_rem - t_call
    t_run <- t_run + t_call

  }

  # Add slurry production rate to output
  dat$slurry_prod_rate <- pars$slurry_prod_rate

  return(dat)

}

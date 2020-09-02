abm_variable <-
function(days, delta_t, pars, warn) {

  # Cannot have no slurry present because is used in all concentration calculations (NTS: could change this)
  pars$slurry_mass[pars$slurry_mass[, 'slurry_mass'] == 0, 'slurry_mass'] <- 1E-6

  # Check for sorted time
  if (is.unsorted(pars$slurry_mass$time)) {
    stop('Column `time` must be sorted when `slurry_mass` is time-dependent, but it not: ',
         pars$slurry_mass$time)
  }

  # Determine slurry removal quantity in each time interval
  # Note final 0--alignment is a bit tricky
  removals <- - c(diff(pars$slurry_mass[, 'slurry_mass']), 0)
  removals[removals < 0] <- 0

  # Remove (combine) consecutive removals because all removals are assumed to be instant
  pars$slurry_mass <- pars$slurry_mass[c(TRUE, !(removals[-1] > 0 & removals[-length(removals)] > 0)), ]

  # And recalculate removals from change in slurry mass
  removals <- - c(diff(pars$slurry_mass[, 'slurry_mass']), 0)
  removals[removals < 0] <- 0

  # Where removal appeared to occur over time, change time to previous value so is instant
  pars$slurry_mass[c(FALSE, which.adj <- diff(pars$slurry_mass$time) > 0 & removals[-length(removals)] > 0), 'time'] <- 
    pars$slurry_mass[c(diff(pars$slurry_mass$time) > 0 & removals[-length(removals)] > 0, FALSE), 'time'] 
  if (any(which.adj) & warn) {
    warning('Non-instant removals given in `slurry_mass` were changed to instant.')
  }

  # Determine variable slurry production rate for each time interval
  slurry_prod_rate_t <- c(diff(pars$slurry_mass[, 'slurry_mass']) / diff(pars$slurry_mass[, 'time']), 0)
  slurry_prod_rate_t[slurry_prod_rate_t < 0] <- 0
  slurry_prod_rate_t[!is.finite(slurry_prod_rate_t)] <- 0

  # Note about time: 1) All simulations start at 0, 2) days must be at least as long as mass data
  add_int <- 0
  t_1 <- NULL
  t_end <- NULL

  if (pars$slurry_mass[1, 'time'] > 0) {
    add_int <- add_int + 1
    t_1 <- 0
    slurry_prod_rate_t <- c(0, slurry_prod_rate_t)
    removals <- c(0, removals)
  }

  if (pars$slurry_mass[1, 'time'] < days) {
    add_int <- add_int + 1
    t_end <- days
    slurry_prod_rate_t <- c(slurry_prod_rate_t, 0)
    removals <- c(removals, 0)
  }

  # Fill in slurry_prod_rate (will be updated for each iteration in loop)
  pars$slurry_prod_rate <- slurry_prod_rate_t[1]

  n_int <- nrow(pars$slurry_mass) + add_int - 1
  t_ints <- diff(c(t_1, pars$slurry_mass[, 'time'], t_end))

  if (any(t_ints < 0)) stop('9x18347hhab')

  # Initial state variable vector
  y <- c(xa = pars$xa_init * pars$slurry_mass[1, 'slurry_mass'], 
         slurry_mass = pars$slurry_mass[1, 'slurry_mass'], 
         Sp = pars$conc_fresh[['Sp']] * pars$slurry_mass[1, 'slurry_mass'], 
         VFA = pars$conc_fresh[['VFA']] * pars$slurry_mass[1, 'slurry_mass'],
         sulfate = SO4_fun(0) * pars$slurry_mass[1, 'slurry_mass'], 
         sulfide = pars$conc_fresh[['S2']] * pars$slurry_mass[1, 'slurry_mass'], 
         CH4_emis_cum = 0, CO2_emis_cum = 0, COD_conv_cum = 0, COD_conv_cum_meth = 0, COD_conv_cum_respir = 0, COD_conv_cum_sr = 0)

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
    t_call <- min(t_ints[i], t_rem)

    # Get number of microbial groups
    n_mic <- length(pars$qhat_opt)

    # Call lsoda and advance if time changes (i.e., if there is not a removal)
    if (t_call > 1E-6) {
      # Need some care with times to make sure t_call is last one in case it is not multiple of delta_t
      times <- sort(unique(round(c(seq(0, t_call, by = min(t_rem, delta_t)), t_call), 5)))

      # Add run time to pars so rates() can use actual time to calculate temp_C and pH
      pars$t_run <- t_run

      # Call up ODE solver
      out <- deSolve::lsoda(y = y, times = times, rates, parms = pars)

      # Extract new state variable vector from last row of lsoda output
      y <- out[nrow(out), 1:(length(pars$qhat_opt) + 11) + 1]

      # Change format of output and drop first (time 0) row (duplicated in last row of previous)
      if (i == 1) {
        out <- data.frame(out)
      } else {
        out <- data.frame(out[-1, , drop = FALSE])
      }

      # Fix some names (NTS: can we sort this out in rates()? I have tried and failed)
      names(out)[1:length(pars$qhat_opt) + 1] <- names(pars$qhat_opt)

      # Change time in output to cumulative time for complete simulation
      out$time <- out$time + t_run

      # Add slurry production rate to output
      out$slurry_prod_rate <- pars$slurry_prod_rate

  
      # Add results to earlier ones
      dat <- rbind(dat, out)
    } else {
      # Empty channel (instantaneous changes at end of day) in preparation for next lsoda call
      # Note: Do not update out here before next iteration
      resid_frac <- 1 - removals[i] / y['slurry_mass']
      resid_xa <- logistic(logit(resid_frac) + pars$resid_enrich)
      
      y[1:n_mic] <- resid_xa * y[1:n_mic]
      y['slurry_mass'] <- resid_frac * y['slurry_mass']
      y['Sp'] <- resid_frac * y['Sp']
      y['VFA'] <- resid_frac * y['VFA']
      y['sulfate'] <- resid_frac * y['sulfate']
      y['sulfide'] <- resid_frac * y['sulfide']
    }
 
    # Update time remaining and total time run so far
    t_rem <- t_rem - t_call
    t_run <- t_run + t_call

    pars$slurry_prod_rate <- slurry_prod_rate_t[i + 1]
   
  }

  return(dat)

}

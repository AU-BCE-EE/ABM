abm_regular <-
function(days, delta_t, pars, starting = NULL, temp_C_fun = temp_C_fun, pH_fun = pH_fun, SO4_fun = SO4_fun) {

  # Cannot have no slurry present because is used in all concentration calculations (NTS: could change this)
  if (pars$slurry_mass == 0) {
    pars$slurry_mass <- 1E-10
  }

  # Sort out emptying interval and related vars
  if (pars$slurry_prod_rate > 0) {
    # Emptying interval (d)
    empty_int <- pars$max_slurry_mass * (1 - pars$resid_frac) / pars$slurry_prod_rate
    # Time of first filling interval (first call may be longer or shorter than rest, storage may be empty or > resid.frac)
    t_int_1 <- (pars$max_slurry_mass - pars$slurry_mass) / pars$slurry_prod_rate
    # Number of filling intervals
    # 1 extra interval for first interval (filling)
    n_int <- ceiling(round((days - t_int_1) / empty_int, 4)) + 1
  } else if (pars$slurry_prod_rate == 0) {
    empty_int <- days
    t_int_1 <- days
    n_int <- 1
  } else stop('slurry_prod_rate problem')

  # Initial state variable vector
  y <- c(xa = pars$xa_init * pars$slurry_mass, 
         slurry_mass = pars$slurry_mass, 
         Sp = pars$conc_fresh[['Sp']] * pars$slurry_mass, 
         VFA = pars$conc_fresh[['VFA']] * pars$slurry_mass, 
         sulfate = SO4_fun(0) * pars$slurry_mass, 
         sulfide = pars$conc_fresh[['S2']] * pars$slurry_mass, 
         CH4_emis_cum = 0, CO2_emis_cum = 0, COD_conv_cum = 0, COD_conv_cum_meth = 0, COD_conv_cum_respir = 0, COD_conv_cum_sr = 0)

  if (!is.null(starting) & is.data.frame(starting)) {
    y['Sp'] <- starting[nrow(starting), 'Sp']
    y['VFA'] <- starting[nrow(starting), 'VFA']
    y['sulfate'] <- starting[nrow(starting), 'sulfate']
    y['sulfide'] <- starting[nrow(starting), 'sulfide']
  }

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
    if (i == 1) {
      t_call <- t_int_1
    } else {
      t_call <- min(empty_int, t_rem)
    }

    # Need some care with times to make sure t_call is last one in case it is not multiple of delta_t
    times <- sort(unique(round(c(seq(0, t_call, by = min(t_rem, delta_t)), t_call), 5)))

    # Add run time to pars so rates() can use actual time to calculate temp_C and pH
    pars$t_run <- t_run

    # Call up ODE solver
    #cat(t_rem, '\n')
    out <- deSolve::lsoda(y = y, times = times, rates, parms = pars, temp_C_fun = temp_C_fun, pH_fun = pH_fun, SO4_fun = SO4_fun)

    # Get number of microbial groups
    n_mic <- length(pars$qhat_opt)

    # Extract new state variable vector from last row
    y <- out[nrow(out), 1:(length(pars$qhat_opt) + 11) + 1]

    # Empty channel (instantaneous changes at end of day) in preparation for next lsoda call
    resid_frac <- pars$resid_frac
    resid_xa <- logistic(logit(resid_frac) + pars$resid_enrich)

    y[1:n_mic] <- resid_xa * y[1:n_mic]
    y['slurry_mass'] <- max(resid_frac * y['slurry_mass'],  1E-10)
    y['Sp'] <- resid_frac * y['Sp']
    y['VFA'] <- resid_frac * y['VFA']
    y['sulfate'] <- resid_frac * y['sulfate']
    y['sulfide'] <- resid_frac * y['sulfide']
    
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

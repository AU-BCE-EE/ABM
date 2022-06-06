abm_regular <- function(days, delta_t, y, pars, temp_fun = temp_fun, pH_fun = pH_fun, SO4_fun = SO4_fun) {

  pars$abm_regular <- TRUE
  
  empty_int <- pars$empty_int

  # Cannot have no slurry present because is used in all concentration calculations (NTS: could change this)
  if (pars$slurry_mass == 0) {
    pars$slurry_mass <- 1E-10
  }

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

  ### NTS: sort out starting in abm.R

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
    out <- deSolve::lsoda(y = y, times = times, rates, parms = pars, temp_fun = temp_fun, pH_fun = pH_fun, SO4_fun = SO4_fun)

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
          empty_names <- grep('^xa\\.|cum|slurry_mass', names(y), value = TRUE, invert = TRUE)
          resid_xa <- logistic(logit(resid_frac) + pars$resid_enrich)
          y[1:n_mic] <- resid_xa * y[1:n_mic]
        } else {
          empty_names <- grep('^xa\\.|cum|slurry_mass|sed$', names(y), value = TRUE, invert = TRUE)
        }
        y[empty_names] <- sapply(y[empty_names], function(x) resid_frac * x) # apply resid_frac to other state variables

      } else {
        warning('Emptying skipped because of low slurry level.')
      }
    }
   
    # Washing, increase slurry mass
    # NTS: needs some work for xa and mix????
    if (wash[i]) {
      y['slurry_mass'] <- y['slurry_mass'] + pars$wash_water

      # And empty again, leaving same residual mass as always
      if (y['slurry_mass'] > pars$resid_mass) {
        resid_wash <- pars$resid_mass / y['slurry_mass']
        y[!grepl('cum', names(y))] <-  y[!grepl('cum', names(y))] * resid_wash
      } else {
        warning('Emptying skipped after adding wash water because of low slurry level.')
      }

      if (pars$rest_d > 0) {
        times <- seq(0, pars$rest_d, delta_t)
        parsr <- pars
        parsr$slurry_prod_rate <- 0
        outr <- deSolve::lsoda(y = y, times = times, rates, parms = parsr, temp_fun = temp_fun, pH_fun = pH_fun, SO4_fun = SO4_fun)
        # Extract new state variable vector from last row
        y <- outr[nrow(outr), 1:(length(y)) + 1]
        # Correct time in outr and combine with main output
        outr[, 'time'] <- outr[, 'time'] + out[nrow(out), 'time']
        out <- rbind(out, outr)
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
    t_rem <- t_rem - t_call - wash[i] * rest_d
    t_run <- t_run + t_call + wash[i] * rest_d

  }

  # Add slurry production rate to output
  dat$slurry_prod_rate <- pars$slurry_prod_rate

  return(dat)

}

abm_regular <- function(days, delta_t, times_regular, y, pars, 
                        temp_C_fun = temp_C_fun, pH_fun = pH_fun, 
                        SO4_inhibition_fun = SO4_inhibition_fun, 
                        conc_fresh_fun = conc_fresh_fun, 
                        xa_fresh_fun = xa_fresh_fun) { 
  
  # Initialize empty interval
  empty_int <- pars$empty_int
  if(empty_int == 0 || is.na(empty_int)) empty_int <- days + 1
  
  if (pars$slurry_mass == 0) pars$slurry_mass <- 1E-10
  
  # Determine wash intervals
  if (!is.na(pars$wash_int) && pars$wash_water > 0) {  
    wash_int <- pars$wash_int
    rest_d <- pars$rest_d
  } else {
    wash_int <- Inf
    rest_d <- 0
  }
  wash_rest_int <- wash_int + rest_d
  
  # Generate intervals
  t_int <- c()
  wash <- c()
  t_nowash <- 0
  while(sum(t_int, wash * rest_d) < days) {
    t_next <- min(wash_int - t_nowash, empty_int, days - sum(t_int, wash * rest_d))
    t_int <- c(t_int, t_next)
    if(t_next == wash_int - t_nowash) {
      wash <- c(wash, TRUE)
      t_nowash <- 0
    } else {
      wash <- c(wash, FALSE)
      t_nowash <- t_nowash + t_next
    }
  }
  
  n_int <- length(t_int)
  
  dat <- NULL
  t_rem <- days
  t_run <- 0
  
  for(i in seq_len(n_int)) {
    
    t_call <- t_int[i]
    times <- sort(unique(round(c(seq(0, t_call, by = min(t_rem, delta_t)), t_call), 5)))
    
    # Local copies to avoid locked binding issues
    y_local <- y
    pars_local <- pars
    pars_local$t_run <- t_run
    pars_local$t_call <- t_call
    pars_local$times <- times
    pars_local$p_idx <- pars_indices(pars_local)
    pars_local$temp_C_fun <- temp_C_fun
    pars_local$pH_fun <- pH_fun
    pars_local$SO4_inhibition_fun <- SO4_inhibition_fun
    pars_local$CTM_cpp <- CTM_cpp
    pars_local$H2SO4_titrate <- H2SO4_titrate
    pars_local$xa_fresh_fun <- xa_fresh_fun
    pars_local$conc_fresh_fun <- conc_fresh_fun
    
    # Call ODE solver
    out <- deSolve::lsoda(y = y_local, times = times, rates_cpp, parms = pars_local)
    
    # Extract last state
    n_state <- length(c(pars$qhat_opt)) + 27
    y <- as.numeric(out[nrow(out), 2:(n_state + 1)])
    
    # Empty store safely
    y <- emptyStore(y, resid_mass = pars_local$resid_mass, resid_enrich = pars_local$resid_enrich)
    y.eff <- y[grepl('_eff$', names(y))]
    y <- y[!grepl('_eff$', names(y))]
    
    # Wash handling
    if(wash[i]) {
      y['slurry_mass'] <- y['slurry_mass'] + pars$wash_water
      y <- emptyStore(y, resid_mass = pars_local$resid_mass, resid_enrich = pars_local$resid_enrich)
      y.eff <- y.eff + y[grepl('_eff$', names(y))]
      y <- y[!grepl('_eff$', names(y))]
      
      if(rest_d > 0) {
        times_rest <- seq(0, rest_d, delta_t)
        pars_rest <- pars_local
        pars_rest$slurry_prod_rate <- 0
        outr <- deSolve::lsoda(y = y, times = times_rest, rates_cpp, parms = pars_rest)
        y <- as.numeric(outr[nrow(outr), 2:(n_state + 1)])
        outr[, 'time'] <- outr[, 'time'] + out[nrow(out), 'time']
        out <- rbind(out, outr)
      }
    }
    
    # Format output
    out <- data.frame(out)
    out[, names(y.eff)] <- 0
    out[nrow(out), names(y.eff)] <- y.eff
    names(out)[1:length(pars$qhat_opt) + 1] <- names(pars$qhat_opt)
    out$time <- out$time + t_run
    
    dat <- dplyr::bind_rows(dat, out)
    t_rem <- t_rem - t_call - wash[i] * rest_d
    t_run <- t_run + t_call + wash[i] * rest_d
  }
  
  if(!is.null(times_regular)) {
    times_regular <- sort(unique(c(times_regular, days)))
    dat <- dat[dat$time %in% times_regular, ]
  }
  
  return(dat)
}

abm_variable <- function(days, delta_t, times = NULL, y, pars, warn = TRUE, 
                         temp_C_fun = temp_C_fun, pH_fun = pH_fun, 
                         SO4_inhibition_fun = SO4_inhibition_fun, 
                         conc_fresh_fun = conc_fresh_fun, 
                         xa_fresh_fun = xa_fresh_fun, 
                         slurry_mass_approx) {
  
  
  
  # Ensure slurry mass never zero
  pars$slurry_mass$slurry_mass[pars$slurry_mass$slurry_mass == 0] <- 1E-10
  
  # Trim and extend slurry_mass for days
  pars$slurry_mass <- pars$slurry_mass[pars$slurry_mass$time <= days, ]
  if (nrow(pars$slurry_mass) == 0 || is.unsorted(pars$slurry_mass$time)) {
    stop("slurry_mass times must be sorted and non-empty")
  }
  
  if (pars$slurry_mass[nrow(pars$slurry_mass), 'time'] < days) {
    last_row <- pars$slurry_mass[nrow(pars$slurry_mass), ]
    last_row$time <- days
    if(ncol(pars$slurry_mass) > 2) last_row[3:ncol(last_row)] <- 0
    pars$slurry_mass <- rbind(pars$slurry_mass, last_row)
  }
  
  # Construct solver times
  if (is.null(times)) times <- seq(0, days, by = delta_t)
  times <- sort(unique(c(times, days, pars$slurry_mass$time)))
  
  # Compute slurry production rates
  slurry_prod_rate_t <- c(0, diff(pars$slurry_mass$slurry_mass) / diff(pars$slurry_mass$time))
  slurry_prod_rate_t[slurry_prod_rate_t < 0 | !is.finite(slurry_prod_rate_t)] <- 0
  
  # Determine removal events
  if (slurry_mass_approx == "late") {
    removals <- -c(0,0,diff(pars$slurry_mass$slurry_mass[-nrow(pars$slurry_mass)])) > 0
  } else if (slurry_mass_approx == "early") {
    removals <- -c(0,diff(pars$slurry_mass$slurry_mass)) > 0
  } else if (slurry_mass_approx == "mid") {
    ir <- which(-c(0,diff(pars$slurry_mass$slurry_mass)) > 0)
    mid_times <- (pars$slurry_mass$time[ir] + pars$slurry_mass$time[ir-1])/2
    mid_slurry <- pars$slurry_mass$slurry_mass[ir]
    pars$slurry_mass <- rbind(pars$slurry_mass, data.frame(time = mid_times, slurry_mass = mid_slurry))
    pars$slurry_mass <- pars$slurry_mass[order(pars$slurry_mass$time), ]
    removals <- -c(0,0,diff(pars$slurry_mass$slurry_mass[-nrow(pars$slurry_mass)])) > 0
    slurry_mass_approx <- "late"
  }
  
  wash_water <- if("wash_water" %in% names(pars$slurry_mass)) pars$slurry_mass$wash_water else rep(0, nrow(pars$slurry_mass))
  
  # Create time intervals
  n_int <- nrow(pars$slurry_mass)
  timelist <- vector("list", n_int)
  cumtime <- numeric(n_int)
  
  for(i in 2:n_int){
    tt <- times[times > pars$slurry_mass$time[i-1] & times <= pars$slurry_mass$time[i]]
    tt <- unique(c(0, tt - max(cumtime[i-1], 0), pars$slurry_mass$time[i] - max(cumtime[i-1],0)))
    timelist[[i]] <- tt
    cumtime[i] <- max(tt) + cumtime[i-1]
  }
  
  dat <- NULL
  t_rem <- days
  t_run <- 0
  
  for(i in 2:n_int){
    t_call <- min(max(timelist[[i]]), t_rem)
    
    # Skip invalid intervals
    if(t_call <= 0) next
    
    # Local copies
    y_local <- y
    pars_local <- pars
    pars_local$t_run <- t_run
    pars_local$slurry_prod_rate <- slurry_prod_rate_t[i]
    pars_local$p_idx <- pars_indices(pars_local)
    pars_local$temp_C_fun <- temp_C_fun
    pars_local$pH_fun <- pH_fun
    pars_local$SO4_inhibition_fun <- SO4_inhibition_fun
    pars_local$CTM_cpp <- CTM_cpp
    pars_local$H2SO4_titrate <- H2SO4_titrate
    pars_local$xa_fresh_fun <- xa_fresh_fun
    pars_local$conc_fresh_fun <- conc_fresh_fun
    
    # Remove slurry if needed
    if(removals[i]){
      j <- if(slurry_mass_approx=="late") i-1 else i
      y_local <- emptyStore(y_local, resid_mass = pars$slurry_mass$slurry_mass[j], resid_enrich = pars$resid_enrich)
      y.eff <- y_local[grepl('_eff$', names(y_local))]
      y_local <- y_local[!grepl('_eff$', names(y_local))]
      
      if(wash_water[i] > 0){
        y_local['slurry_mass'] <- y_local['slurry_mass'] + wash_water[i]
        y_local <- emptyStore(y_local, resid_mass = pars$slurry_mass$slurry_mass[j], resid_enrich = pars$resid_enrich)
        y.eff <- y.eff + y_local[grepl('_eff$', names(y_local))]
        y_local <- y_local[!grepl('_eff$', names(y_local))]
      }
    }
    
    # Solve ODE
    out <- deSolve::lsoda(y = y_local, times = timelist[[i]], rates_cpp, parms = pars_local)
    if(i>2) out <- out[-1, , drop=FALSE]
    out <- data.frame(out)
    
    names(out)[2:(length(pars$qhat_opt)+1)] <- names(pars$qhat_opt)
    y <- unlist(out[nrow(out), 2:(length(pars$qhat_opt)+1)])
    
    out$time <- out$time + t_run
    
    # Add effluent
    yy <- emptyStore(y, resid_mass = 0, resid_enrich = 0)
    y.eff <- 0*yy[grepl('_eff$', names(yy))]
    out[, names(y.eff)] <- 0
    out[nrow(out), names(y.eff)] <- y.eff
    
    dat <- dplyr::bind_rows(dat, out)
    t_rem <- t_rem - t_call
    t_run <- t_run + t_call
  }
  
  # Drop times that were added artificially
  times_orig <- sort(unique(c(times, days)))
  droptimes <- times[!times %in% times_orig]
  dat <- dat[!dat$time %in% droptimes, ]
  
  return(dat)
}

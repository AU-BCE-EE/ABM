
abm <- function(
  days = 365,                                # Number of days to run
  delta_t = 1,                               # Time step for output
  times = NULL,
  mng_pars = list(slurry_prod_rate = 5700,   # kg/d
                  slurry_mass = 39000,       # Initial slurry mass (kg) 
                  storage_depth = 0.6,       # Storge structure depth, assued to be maximum slurry depth (m)
                  resid_depth = 0.05,        # Residual slurry depth after emptying (m)
                  area = 715,                # Tank/pit slurry surface area (assume vertical sides, includes slurry underneath walking path) (m2)
                  empty_int = 42,            # Fixed emptying interval (days)
                  temp_C = 20,
                  wash_water = 0,            
                  wash_int = NA,
                  rest_d = 0,
                  cover = 'none',
                  resid_enrich = 1,
                  scale = c(ks_coefficient = 1, qhat_opt = 1, xa_fresh = 1, yield = 1, alpha_opt = 1)),
  man_pars = man_parsx,
  init_pars = list(conc_init = man_pars$conc_fresh),
  grp_pars = grp_parsx,
  mic_pars = mic_parsx,
  chem_pars = chem_parsx,
  ctrl_pars = list(respir = TRUE,
                   pH_inhib = FALSE, 
                   approx_method = c(temp = 'linear', pH = 'linear', slurry_mass = 'early'), 
                   par_key = '\\.',
                   rates_calc = 'instant'),
  add_pars = NULL,
  pars = NULL,
  startup = 0,                                # Number of times complete simulation should be run before returning results
  starting = NULL,                            # Output from previous simulation to be starting condition for new one
  value = 'ts',                               # Type of output
  warn = TRUE) {

  # If startup repetitions are requested, repeat some number of times before returning results
  # This is done in a for loop in abm_repeat()
  #if (startup > 0) {
  #  cat('\nStartup run ')
  #  out <- do.call(abm_repeat, as.list(environment()))
  #  return(out)
  #} 

  # Sort out parameters, ultimately packaging all parameters into a single list pars
  pars <- packPars(mng_pars = mng_pars,
                   man_pars = man_pars,
                   init_pars = init_pars,
                   grp_pars = grp_pars,
                   mic_pars = mic_pars,
                   chem_pars = chem_pars,
                   ctrl_pars = ctrl_pars,
                   add_pars = add_pars,
                   pars = pars,
                   starting = starting,
                   days = days)

  # Create time-variable functions
  # The element pars$x must either be a numeric constant for constant temperature or pH, or a data frame with time (column 1) and variable value (column 2)
  temp_C_fun <- makeTimeFunc(pars$temp_C, approx_method = pars$approx_method['temp'])
  pH_fun <- makeTimeFunc(pars$pH, approx_method = pars$approx_method['pH'])

  ## Inhibition function
  #SO4_inhib_fun <- approxfun(c(0, 0.09, 0.392251223, 0.686439641, 1.046003263, 1.470942088, 1.961256117, 4), 
  #                           c(1, 1, 0.85, 0.3172, 0.29, 0.1192, 0.05, 0.001), rule = 2)

  # Create initial state variable vector
  y <- makeInitState(pars, starting) 

  if (any(pars$conc_fresh[['VSd']] > 2e-10) & any(pars$conc_fresh[names(pars$conc_fresh) %in% names(pars$A)[!names(pars$A) %in% c('VSd','urea')]] > 2)) {
    stop('Cannot have both VSd and other organic matter components being above 0')
  }
  
  # hard wired parameters, that will not change during a rates call
  # and was moved from rates to here to speed up model. These parameters are added in the hard_pars().
  pars <- hard_pars(pars)

  if (is.numeric(pars$slurry_mass)) {
    # Option 1: Fixed slurry production rate, regular emptying schedule
    dat <- abm_regular(days = days, 
                       delta_t = delta_t, 
                       times_regular = times, 
                       y = y, 
                       pars = pars, 
                       starting = starting, 
                       temp_C_fun = temp_C_fun, 
                       pH_fun = pH_fun)
  } else if (is.data.frame(pars$slurry_mass)) {
    # Option 2: Everything based on given slurry mass vs. time
    dat <- abm_variable(days = days, 
                        delta_t = delta_t, 
                        times = times, 
                        y = y, 
                        pars = pars, 
                        warn = warn, 
                        temp_C_fun = temp_C_fun, 
                        pH_fun = pH_fun, 
                        slurry_mass_approx = pars$approx_method['slurry_mass'])
  } 

  #colnames(dat) <- gsub("conc_fresh.","conc_fresh_", colnames(dat))
  
  ## Calculate rates etc. for output, from state variables
  ## NTS: check units, use dens???
  #if(rates_calc != 'instant'){
  #  dat$rNH3 <- c(0, diff(dat$NH3_emis_cum))/c(1, diff(dat$time))
  #  dat$NH3_emis_rate <- c(0, diff(dat$NH3_emis_cum))/c(1, diff(dat$time))

  #  dat$rN2O <- c(0, diff(dat$N2O_emis_cum))/c(1, diff(dat$time))
  #  dat$N2O_emis_rate <- c(0, diff(dat$N2O_emis_cum))/c(1, diff(dat$time))
  #  
  #  dat$rCH4 <- c(0, diff(dat$CH4_emis_cum))/c(1, diff(dat$time))
  #  dat$CH4_emis_rate <- c(0, diff(dat$CH4_emis_cum))/c(1, diff(dat$time))
 
  #  dat$rCO2 <- c(0, diff(dat$CO2_emis_cum))/c(1, diff(dat$time))
  #  dat$CO2_emis_rate <- c(0, diff(dat$CO2_emis_cum))/c(1, diff(dat$time))
  #}
  
  # NTS: Put below stuff and more into cleanOut() function
  # Replace . in names with _
  names(dat) <- gsub('\\.', '_', names(dat))
  # Add slurry depth
  dat$slurry_depth <- dat$slurry_mass / pars$area / pars$dens

  # Return results
  # Average only
  if (substring(value, 1, 3) == 'sum') {
    return(summ)
  }
  
  # ts = time series
  if (value == 'ts') {
    return(dat)
  }
  
  if (substring(value, 1, 3) == 'eff') {
    return(eff)
  }

  # NTS: can include in everything below too
  # Or everything
  return(list(pars = pars, ts = dat, summ = summ))

}


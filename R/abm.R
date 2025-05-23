
abm <- function(
  days = 365,                                # Number of days to run
  delta_t = 1,                               # Time step for output
  times = NULL,                              # Optional vector of times for output
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

  # Sort out parameters, packaging all parameters into a single list pars, adding some others
  pars <- packPars(mng_pars = mng_pars,
                   man_pars = man_pars,
                   init_pars = init_pars,
                   grp_pars = grp_pars,
                   mic_pars = mic_pars,
                   chem_pars = chem_pars,
                   ctrl_pars = ctrl_pars,
                   add_pars = add_pars,
                   pars = pars,
                   starting = starting)
  
  # If startup repetitions are requested, repeat some number of times before returning results
  # This is done in a for loop in abm_repeat()
  if (startup > 0) {
    cat('\nStartup run ')
    out <- abmStartup(days = days,
                      delta_t = delta_t,
                      times = times,
                      pars = pars,
                      startup = startup,
                      starting = starting,
                      value = value,
                      warn = warn)
    return(out)
  } 


  # Create time-variable functions
  # The element pars$x must either be a numeric constant for constant temperature or pH, or a data frame with time (column 1) and variable value (column 2)
  temp_C_fun <- makeTimeFunc(pars$temp_C, approx_method = pars$approx_method['temp'])
  pH_fun <- makeTimeFunc(pars$pH, approx_method = pars$approx_method['pH'])

  # Create initial state variable vector
  y <- makeInitState(pars) 

  if (is.numeric(pars$slurry_mass)) {
    # Option 1: Fixed slurry production rate, regular emptying schedule
    dat <- abm_regular(days = days, 
                       delta_t = delta_t, 
                       times_regular = times, 
                       y = y, 
                       pars = pars, 
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

  # Clean up and possibly extend output
  dat <- cleanOutput(dat, pars, addcols = TRUE, addconcs = TRUE)

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

  # Or everything
  return(list(pars = pars, ts = dat, summ = summ))

}



abm <- function(
  days = 365,                                # Number of days to run
  delta_t = 1,                               # Time step for output
  times = NULL,                              # Optional vector of times for output
  mng_pars = NULL,
  man_pars = NULL,
  init_pars = list(conc_init =  c(man_pars$comp_fresh, man_pars$VFA_fresh)),
  grp_pars = NULL,
  mic_pars = NULL,
  sub_pars = NULL,
  chem_pars = NULL,
  inhib_pars = NULL,
  mt_pars = NULL,
  ctrl_pars = list(respir = TRUE,
                   pH_inhib = FALSE, 
                   approx_method = 'early',
                   par_key = '\\.',
                   rates_calc = 'instant'),
  var_pars = list(var = NULL),
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
                   sub_pars = sub_pars,
                   chem_pars = chem_pars,
                   inhib_pars = inhib_pars,
                   ctrl_pars = ctrl_pars,
                   var_pars = var_pars,
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


  # Create initial state variable vector
  y <- makeInitState(pars) 

  if (is.null(pars$var)) {
    # Option 1: Fixed slurry production rate, regular emptying schedule
    dat <- abmReg(days = days, 
                       delta_t = delta_t, 
                       times_regular = times, 
                       y = y, 
                       pars = pars)
  } else if (inherits(pars$var, 'data.frame')) {
    # Option 2: Everything based on given slurry mass vs. time
    dat <- abmVar(days = days, 
                        delta_t = delta_t, 
                        times = times, 
                        y = y, 
                        pars = pars, 
                        warn = warn)
  } else {
    stop('pars_var must be NULL or a data frame')
  }

  # Clean up and possibly extend output
  dat <- cleanOutput(dat, pars, addcols = TRUE, addconcs = TRUE, cumeff = TRUE)

  # Check COD balance
  codbal <- checkCOD(dat = dat, grps = pars$grps, subs = pars$subs, COD_conv = pars$COD_conv, rtol = 0.01)

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


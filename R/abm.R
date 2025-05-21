#abm
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
                  wash_water = 75000,            
                  wash_int = NA,
                  rest_d = 5,
                  cover = 'none',
                  resid_enrich = 0.9,
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

  colnames(dat) <- gsub("conc_fresh.","conc_fresh_", colnames(dat))
  
  # Calculate concentrations where relevant
  #conc.names <- names(dat)[!grepl('conc|time|slurry_mass|inhib|qhat|CH4_emis_cum', names(dat))]
  mic_names <- pars$grps
  eff_names <- names(dat[grepl("_eff$", names(dat))])
  eff_conc_names <- eff_names[eff_names != "slurry_mass_eff"]
  conc_names <-  c('TAN', 'xa_aer', 'xa_bac', 'xa_dead', 'urea', 'RFd', 'iNDF', 'ash', 'VSd', 'starch', 'Cfat', 'CPs', 'CPf', 'VFA', 'sulfide', 'sulfate', mic_names)
  dat_conc <- dat[, conc_names]/(dat$slurry_mass)
  dat_eff_conc <- dat[, eff_conc_names]/(dat$slurry_mass_eff)
  names(dat_conc) <- paste0(names(dat_conc), '_conc')
  names(dat_eff_conc) <- paste0(names(dat_eff_conc), '_conc')
  dat <- cbind(dat, dat_conc, dat_eff_conc)

  # Add temperature and pH
  dat$temp_C <- temp_C_fun(dat$time)
  if (is.numeric(pars$pH) | is.data.frame(pars$pH)) {
    dat$pH <- pH_fun(dat$time)
  } else if (pars$pH == 'calc'){
    dat$pH <- H2SO4_titrat(dat$sulfate_conc, class_anim = "pig")$pH
  } else {
    stop('Problem with pH input (bee721)')
  }
  
  # Calculate rates etc. for output, from state variables
  # NTS: check units, use dens???

  dat$NH3_emis_rate <- dat$NH3_emis_rate_pit + dat$NH3_emis_rate_floor
  
  if(rates_calc != 'instant'){
  dat$rNH3 <- c(0, diff(dat$NH3_emis_cum))/c(1, diff(dat$time))
  dat$NH3_emis_rate <- c(0, diff(dat$NH3_emis_cum))/c(1, diff(dat$time))

  dat$rN2O <- c(0, diff(dat$N2O_emis_cum))/c(1, diff(dat$time))
  dat$N2O_emis_rate <- c(0, diff(dat$N2O_emis_cum))/c(1, diff(dat$time))
  
  dat$rCH4 <- c(0, diff(dat$CH4_emis_cum))/c(1, diff(dat$time))
  dat$CH4_emis_rate <- c(0, diff(dat$CH4_emis_cum))/c(1, diff(dat$time))
 
  dat$rCO2 <- c(0, diff(dat$CO2_emis_cum))/c(1, diff(dat$time))
  dat$CO2_emis_rate <- c(0, diff(dat$CO2_emis_cum))/c(1, diff(dat$time))
  }
  
  dat$NH3_emis_rate_slurry <- dat$NH3_emis_rate / (dat$slurry_mass / 1000)
  dat$NH3_flux <- dat$NH3_emis_rate / pars$area
  dat$N2O_emis_rate_slurry <- dat$N2O_emis_rate / (dat$slurry_mass / 1000)
  dat$N2O_flux <- dat$N2O_emis_rate / pars$area
  dat$CH4_emis_rate_slurry <- dat$CH4_emis_rate / (dat$slurry_mass / 1000)
  dat$CH4_flux <- dat$CH4_emis_rate / pars$area
  dat$CO2_emis_rate_slurry <- dat$CO2_emis_rate / (dat$slurry_mass / 1000)
  dat$CO2_flux <- dat$CO2_emis_rate / pars$area
  
  ## NTS: Add others, e.g., mu
  # Calculate COD/VS flows
  # First concentrations in g/kg

  dat$dCOD_conc_fresh <- dat$conc_fresh_VFA + dat$conc_fresh_xa_aer + dat$conc_fresh_xa_bac + dat$conc_fresh_xa_dead + dat$conc_fresh_RFd + dat$conc_fresh_starch + dat$conc_fresh_CPs + dat$conc_fresh_CPf + dat$conc_fresh_Cfat + dat$conc_fresh_VSd + rowSums(dat[, grepl("xa_fresh_", colnames(dat))])
  dat$COD_conc_fresh <- dat$conc_fresh_VFA + dat$conc_fresh_xa_aer + dat$conc_fresh_xa_bac + dat$conc_fresh_xa_dead + dat$conc_fresh_RFd + dat$conc_fresh_starch + dat$conc_fresh_CPs + dat$conc_fresh_CPf + dat$conc_fresh_Cfat + dat$conc_fresh_VSd + rowSums(dat[, grepl("xa_fresh_", colnames(dat))]) + dat$conc_fresh_iNDF
  dat$ndCOD_conc_fresh <- dat$COD_conc_fresh - dat$dCOD_conc_fresh
  dat$C_conc_fresh <- dat$conc_fresh_VFA / pars$COD_conv[['C_VFA']] + dat$conc_fresh_xa_aer / pars$COD_conv[['C_xa']] + dat$conc_fresh_xa_bac / pars$COD_conv[['C_xa']] + dat$conc_fresh_xa_dead / pars$COD_conv[['C_xa']] + dat$conc_fresh_RFd / pars$COD_conv[["C_RFd"]] +
                      dat$conc_fresh_starch / pars$COD_conv[['C_starch']] + dat$conc_fresh_CPs / pars$COD_conv[['C_CP']] + dat$conc_fresh_CPf / pars$COD_conv[['C_CP']] + dat$conc_fresh_Cfat / pars$COD_conv[['C_Cfat']] +
                      dat$conc_fresh_VSd / pars$COD_conv[['C_VSd']] + dat$conc_fresh_iNDF / pars$COD_conv[['C_iNDF']] + rowSums(dat[, grepl("xa_fresh_", colnames(dat))]) / pars$COD_conv[['C_xa']] +
                      dat$conc_fresh_urea / pars$COD_conv[['C_N_urea']]
  
  # ndCOD is almost conserved, same everywhere always - some problem here with water evap and water precipitation?
  dat$ndCOD_conc <- ndCOD_conc <- dat$COD_conc_fresh - dat$dCOD_conc_fresh 
  dat$dCOD_conc <- dCOD_conc <- dat$xa_aer_conc + dat$xa_bac_conc + dat$xa_dead_conc + dat$RFd_conc + dat$Cfat_conc + dat$CPs_conc + dat$CPf_conc + dat$starch_conc + dat$VFA_conc + dat$VSd_conc + rowSums(dat[, paste0(mic_names, '_', 'conc'), drop = FALSE])
  dat$COD_conc <- COD_conc <- ndCOD_conc + dCOD_conc
  
  dat$VS_conc <- dat$COD_conc/pars$COD_conv[['VS']]

  dat$C_conc <- dat$VFA_conc / pars$COD_conv[['C_VFA']] + dat$xa_aer_conc / pars$COD_conv[['C_xa']] + dat$xa_bac_conc / pars$COD_conv[['C_xa']] + dat$xa_dead_conc / pars$COD_conv[['C_xa']] + dat$RFd_conc / pars$COD_conv[["C_RFd"]] +
                dat$starch_conc / pars$COD_conv[['C_starch']] + dat$CPs_conc / pars$COD_conv[['C_CP']] + dat$CPf_conc / pars$COD_conv[['C_CP']] + dat$Cfat_conc / pars$COD_conv[['C_Cfat']] +
                dat$VSd_conc / pars$COD_conv[['C_VSd']] + dat$iNDF_conc / pars$COD_conv[['C_iNDF']] + rowSums(dat[, paste0(mic_names, '_', 'conc'), drop = FALSE]) / pars$COD_conv[['C_xa']] +
                dat$urea_conc / pars$COD_conv[['C_N_urea']]
  
  dat$C_eff_conc <- dat$VFA_eff_conc / pars$COD_conv[['C_VFA']] + dat$xa_aer_eff_conc / pars$COD_conv[['C_xa']] + dat$xa_bac_eff_conc / pars$COD_conv[['C_xa']] + dat$xa_dead_eff_conc / pars$COD_conv[['C_xa']] + dat$RFd_eff_conc / pars$COD_conv[["C_RFd"]] +
                dat$starch_eff_conc / pars$COD_conv[['C_starch']] + dat$CPs_eff_conc / pars$COD_conv[['C_CP']] + dat$CPf_eff_conc / pars$COD_conv[['C_CP']] + dat$Cfat_eff_conc / pars$COD_conv[['C_Cfat']] +
                dat$VSd_eff_conc / pars$COD_conv[['C_VSd']] + dat$iNDF_eff_conc / pars$COD_conv[['C_iNDF']] + rowSums(dat[, paste0(mic_names, '_', 'conc'), drop = FALSE]) / pars$COD_conv[['C_xa']] +
                dat$urea_eff_conc / pars$COD_conv[['C_N_urea']]
  
  # still miss to incorp N from biomass here, 
  # specifically the N stored in xa_bac needs to be transfered to the CP pool when it degrades.  
  
  dat$Ninorg_conc_fresh <- dat$conc_fresh_urea + dat$conc_fresh_TAN
  dat$Norg_conc_fresh <- dat$conc_fresh_CPs / pars$COD_conv[['CP_N']] + dat$conc_fresh_CPf / pars$COD_conv[['CP_N']] + dat$conc_fresh_xa_bac / pars$COD_conv[['N_xa']] + dat$conc_fresh_xa_aer / pars$COD_conv[['N_xa']]# Note that the input CP is typically measured and therefore includes.
  dat$N_conc_fresh <- dat$Ninorg_conc_fresh + dat$Norg_conc_fresh
  dat$Ninorg_conc <- Ninorg_conc <- dat$urea_conc + dat$TAN_conc
  dat$Norg_conc <- Norg_conc <- dat$CPs_conc / pars$COD_conv[['CP_N']] + dat$CPf_conc / pars$COD_conv[['CP_N']] + dat$xa_bac_conc / pars$COD_conv[['N_xa']] + dat$xa_aer_conc / pars$COD_conv[['N_xa']]
  dat$N_conc <- N_conc <- dat$Ninorg_conc + dat$Norg_conc
  dat$Ninorg_eff_conc <- Ninorg_eff_conc <- dat$urea_eff_conc + dat$TAN_eff_conc
  dat$N_eff_conc <- Ninorg_eff_conc + dat$CPs_eff_conc / pars$COD_conv[['CP_N']] + dat$CPf_eff_conc / pars$COD_conv[['CP_N']] +dat$xa_bac_eff_conc / pars$COD_conv[['N_xa']] + dat$xa_aer_eff_conc / pars$COD_conv[['N_xa']]
  
  # And flows in g/d
  dat$dCOD_load_rate <- dat$dCOD_conc_fresh * dat$slurry_prod_rate
  dat$ndCOD_load_rate <- dat$ndCOD_conc_fresh * dat$slurry_prod_rate
  dat$VS_load_rate <- dat$COD_load_rate / pars$COD_conv[['VS']]
  #dat$C_load_rate <- dat$C_conc_fresh * dat$slurry_prod_rate
  dat$slurry_load_rate <- dat$slurry_prod_rate/1000 # m3
  
  dat$Ninorg_load_rate <- dat$Ninorg_conc_fresh * dat$slurry_prod_rate
  dat$Norg_load_rate <- dat$Norg_conc_fresh * dat$slurry_prod_rate
  #dat$N_load_rate <- dat$N_conc_fresh * dat$slurry_prod_rate

  # NTS: These need to be deleted after replacing desired ones with output from lsoda() call
  ## Cumulative flow in g
  ## NTS still need to have an instant COD_load_rate here, we just need to add it in output of rates to calculate the "true" COD_load_cum 
  #dat$COD_load_cum <- cumsum(dat$COD_load_rate * c(0, diff(dat$time))) + dat$COD_conc[1] * dat$slurry_mass[1]
  #dat$dCOD_load_cum <- cumsum(dat$dCOD_load_rate * c(0, diff(dat$time))) + dat$dCOD_conc[1]* dat$slurry_mass[1]
  #dat$ndCOD_load_cum <- cumsum(dat$ndCOD_load_rate * c(0, diff(dat$time))) + dat$ndCOD_conc[1]* dat$slurry_mass[1]
  #dat$VS_load_cum <- cumsum(dat$VS_load_rate * c(0, diff(dat$time))) + dat$VS_conc[1]* dat$slurry_mass[1]
  #dat$C_load_cum <- cumsum(dat$C_load_rate * c(0, diff(dat$time))) + dat$C_conc[1] * dat$slurry_mass[1]
  #dat$slurry_load_cum <- (cumsum(dat$slurry_load_rate * c(0, diff(dat$time))) + dat$slurry_mass[1]/1000)
  #
  #
  #dat$Ninorg_load_cum <- cumsum(dat$Ninorg_load_rate * c(0, diff(dat$time))) + dat$Ninorg_conc[1] * dat$slurry_mass[1]
  #dat$Norg_load_cum <- cumsum(dat$Norg_load_rate * c(0, diff(dat$time))) + dat$Norg_conc[1] * dat$slurry_mass[1]
  #dat$N_load_cum <- cumsum(dat$N_load_rate * c(0, diff(dat$time))) + dat$N_conc[1] * dat$slurry_mass[1]

  # And relative emission
  # g CH4/g COD in
  
  dat$CH4_emis_rate_COD <- dat$CH4_emis_rate / dat$COD_load_rate
  #dat$CH4_emis_rate_dCOD <- dat$CH4_emis_rate / dat$dCOD_load_rate
  #dat$CH4_emis_rate_VS <- dat$CH4_emis_rate / dat$VS_load_rate
  dat$CH4_emis_rate_C <- dat$CH4_emis_rate / dat$C_load_rate
  dat$CH4_emis_cum_COD <- dat$CH4_emis_cum / dat$COD_load_cum
  #dat$CH4_emis_cum_dCOD <- dat$CH4_emis_cum / dat$dCOD_load_cum
  #dat$CH4_emis_cum_VS <- dat$CH4_emis_cum / dat$VS_load_cum
  dat$CH4_emis_cum_C <- dat$CH4_emis_cum / dat$C_load_cum
  #dat$CH4_emis_cum_slurry <- dat$CH4_emis_cum / dat$slurry_load_cum # g/ m3
  
  # g NH3 N/g N
  dat$NH3_emis_rate_Ninorg <- dat$NH3_emis_rate / dat$Ninorg_load_rate
  dat$NH3_emis_rate_Norg <- dat$NH3_emis_rate / dat$Norg_load_rate
  dat$NH3_emis_rate_N <- dat$NH3_emis_rate / dat$N_load_rate
  #dat$NH3_emis_cum_Ninorg <- dat$NH3_emis_cum / dat$Ninorg_load_cum
  #dat$NH3_emis_cum_Norg <- dat$NH3_emis_cum / dat$Norg_load_cum
  dat$NH3_emis_cum_N <- dat$NH3_emis_cum / dat$N_load_cum
  
  # g N2O N/g N
  dat$N2O_emis_rate_Ninorg <- dat$N2O_emis_rate / dat$Ninorg_load_rate
  dat$N2O_emis_rate_Norg <- dat$N2O_emis_rate / dat$Norg_load_rate
  dat$N2O_emis_rate_N <- dat$N2O_emis_rate / dat$N_load_rate
  #dat$N2O_emis_cum_Ninorg <- dat$N2O_emis_cum / dat$Ninorg_load_cum
  #dat$N2O_emis_cum_Norg <- dat$N2O_emis_cum / dat$Norg_load_cum
  dat$N2O_emis_cum_N <- dat$N2O_emis_cum / dat$N_load_cum
  
  # Same for CO2
  dat$CO2_emis_rate_COD <- dat$CO2_emis_rate / dat$COD_load_rate
  dat$CO2_emis_rate_dCOD <- dat$CO2_emis_rate / dat$dCOD_load_rate
  dat$CO2_emis_rate_VS <- dat$CO2_emis_rate / dat$VS_load_rate
  dat$CO2_emis_rate_C <- dat$CO2_emis_rate / dat$C_load_rate
  dat$CO2_emis_cum_COD <- dat$CO2_emis_cum / dat$COD_load_cum
  #dat$CO2_emis_cum_dCOD <- dat$CO2_emis_cum / dat$dCOD_load_cum
  #dat$CO2_emis_cum_VS <- dat$CO2_emis_cum / dat$VS_load_cum
  dat$CO2_emis_cum_C <- dat$CO2_emis_cum / dat$C_load_cum
  #dat$CO2_emis_cum_slurry <- dat$CO2_emis_cum / dat$slurry_load_cum
  
  # Apparent COD conversion fraction
  #dat$f_COD_CH4_rate <-  dat$CH4_emis_rate * pars$COD_conv[['CH4']] / dat$COD_load_rate
  dat$f_COD_CH4_cum <-  dat$COD_conv_cum_meth / dat$COD_load_cum
  dat$f_COD_respir_cum <-  dat$COD_conv_cum_respir / dat$COD_load_cum
  dat$f_COD_sr_cum <-  dat$COD_conv_cum_sr / dat$COD_load_cum
  
  # Replace . in names with _
  names(dat) <- gsub('\\.', '_', names(dat))
  
  Ninorg_load <- dat$Ninorg_load_cum[nrow(dat)]
  Norg_load <- dat$Norg_load_cum[nrow(dat)]
  N_load <- dat$N_load_cum[nrow(dat)]
  
  COD_load <- dat$COD_load_cum[nrow(dat)]
  dCOD_load <- dat$dCOD_load_cum[nrow(dat)]
  ndCOD_load <- dat$ndCOD_load_cum[nrow(dat)]
  VS_load <- dat$VS_load_cum[nrow(dat)]
  C_load <- dat$C_load_cum[nrow(dat)]
  slurry_load <- dat$slurry_load_cum[nrow(dat)]
  
  NH3_emis_cum <- dat$NH3_emis_cum[nrow(dat)]
  NH3_emis_rate <- NH3_emis_cum / (dat$time[nrow(dat)] - dat$time[1])
  NH3_emis_Ninorg <- NH3_emis_cum / Ninorg_load
  NH3_emis_Norg <- NH3_emis_cum / Norg_load
  NH3_emis_N <- NH3_emis_cum / N_load
  
  N2O_emis_cum <- dat$N2O_emis_cum[nrow(dat)]
  N2O_emis_rate <- N2O_emis_cum / (dat$time[nrow(dat)] - dat$time[1])
  N2O_emis_Ninorg <- N2O_emis_cum / Ninorg_load
  N2O_emis_Norg <- N2O_emis_cum / Norg_load
  N2O_emis_N <- N2O_emis_cum / N_load
  
  CH4_emis_cum <- dat$CH4_emis_cum[nrow(dat)]
  CH4_emis_rate <- CH4_emis_cum / (dat$time[nrow(dat)] - dat$time[1])
  CH4_emis_COD <- CH4_emis_cum / COD_load
  CH4_emis_dCOD <- CH4_emis_cum / dCOD_load
  CH4_emis_VS <- CH4_emis_cum / VS_load
  CH4_emis_C <- CH4_emis_cum / C_load
  CH4_emis_slurry <- CH4_emis_cum / slurry_load
  
  CO2_emis_cum <- dat$CO2_emis_cum[nrow(dat)] - dat$CO2_emis_cum[1]
  CO2_emis_rate <- CO2_emis_cum / (dat$time[nrow(dat)] - dat$time[1])
  CO2_emis_COD <- CO2_emis_cum / COD_load
  CO2_emis_dCOD <- CO2_emis_cum / dCOD_load
  CO2_emis_VS <- CO2_emis_cum / VS_load
  CO2_emis_C <- CO2_emis_cum / C_load
  CO2_emis_slurry <- CO2_emis_cum / slurry_load
  
  COD_conv_meth <- dat$COD_conv_cum_meth[nrow(dat)] - dat$COD_conv_cum_meth[1]
  COD_conv_respir <- dat$COD_conv_cum_respir[nrow(dat)] - dat$COD_conv_cum_respir[1]
  COD_conv_sr <- dat$COD_conv_cum_sr[nrow(dat)] - dat$COD_conv_cum_sr[1]
  f_COD_CH4 <- COD_conv_meth / COD_load
  f_COD_respir <- COD_conv_respir / COD_load
  f_COD_sr <- COD_conv_sr / COD_load
  
  summ <- c(C_load = C_load, COD_load = COD_load, dCOD_load = dCOD_load, ndCOD_load = ndCOD_load, VS_load = VS_load, Ninorg_load = Ninorg_load, Norg_load = Norg_load, 
            N_load = N_load, NH3_emis_cum = NH3_emis_cum, NH3_emis_rate = NH3_emis_rate, NH3_emis_Ninorg = NH3_emis_Ninorg, NH3_emis_Norg = NH3_emis_Norg,
            NH3_emis_N = NH3_emis_N, N2O_emis_cum = N2O_emis_cum, N2O_emis_rate = N2O_emis_rate, N2O_emis_Ninorg = N2O_emis_Ninorg, N2O_emis_Norg = N2O_emis_Norg,
            N2O_emis_N = N2O_emis_N, CH4_emis_cum = CH4_emis_cum, CH4_emis_rate = CH4_emis_rate,
            CH4_emis_COD = CH4_emis_COD, CH4_emis_dCOD = CH4_emis_dCOD, CH4_emis_VS = CH4_emis_VS, CH4_emis_C = CH4_emis_C, CH4_emis_slurry = CH4_emis_slurry,
            CO2_emis_cum = CO2_emis_cum, CO2_emis_rate = CO2_emis_rate, CO2_emis_COD = CO2_emis_COD, CO2_emis_dCOD = CO2_emis_dCOD, CO2_emis_VS = CO2_emis_VS, CO2_emis_C = CO2_emis_C, CO2_emis_slurry = CO2_emis_slurry,
            COD_conv_meth = COD_conv_meth, COD_conv_respir = COD_conv_respir, COD_conv_sr = COD_conv_sr,
            f_COD_CH4 = f_COD_CH4, f_COD_respir = f_COD_respir, f_COD_sr = f_COD_sr) 

  # Add slurry depth
  dat$slurry_depth <- dat$slurry_mass / pars$area / pars$dens

  # Get effluent-only lines
  eff <- dat[dat$slurry_mass_eff > 0, ]
  eff$slurry_mass_eff_cum <- cumsum(eff$slurry_mass_eff)

  # Check slurry mass, warn if needed
  if (any(dat$slurry_mass > pars$max_slurry_mass)) {
    warning('Maximum slurry mass exceeded.\nCheck output.')
  }
 
  # print message with inputs
  if(is.data.frame(pars$slurry_mass)){
    message('arguments overwritten by slurry_mass: slurry_prod_rate, empty_int, resid_depth, wash_water, wash_int, rest_d')
  }
  
  message(paste0('rain = ', pars$rain,' kg/m2/day'))
  message(paste0('evaporation = ', round(pars$evap,2),' kg/m2/day'))
  
  if (rates_calc != 'instant' & !is.null(times)) {
    warning(paste0('rates_calc is ', rates_calc,' but specifc output times are given. 
                   It is recommended to change rates_calc to instant'))
  }
  # Return results
  # Average only
  if (substring(value, 1, 3) == 'sum') return(summ)
  # ts = time series
  if (value == 'ts') return(dat)
  if (substring(value, 1, 3) == 'eff') return(eff)

  # NTS: create and return cumulative effluent here?
  # NTS: remove duplicate times from emptying (actually only emptying times should show up in effluent right?)
  # NTS: can include in everything below too
  # Or everything
  
  return(list(pars = pars, ts = dat, summ = summ))

}


atm <-
function(
  days = 365,                               # Number of days to run
  delta_t = 1,                              # Time step for output
  mng_pars = list(slurry_prod_rate = 1000,  # kg/d
                 slurry_mass = 0,           # Initial slurry mass (kg)
                 max_slurry_mass = 33333,   # Maximum slurry mass (kg), default makes for 30 d cycle
                 resid_frac = 0.10,         # Residual slurry fraction after emptying
                 area = 11,                 # Based on 3 m depth (m2)
                 temp_C = 23),
  man_pars = list(conc_fresh = list(S2 = 0.0, SO4 = 0.2, TAN = 1.0, VFA = 4.2, Sp = 65, COD = 160),
                 pH = 7), # Note list not c() so SO4 can be data frame
  grp_pars = list(yield = c(default = 0.04, sr1 = 0.065),
                  xa_fresh = c(default = 0.001, sr1 = 0.001),
                  xa_init = c(all = 0.01),
                  decay_rate = c(all = 0.02),
                  ks_coefficient = c(m1 = 0.5, m2 = 1.5, m3 = 1.0, p1 = 1.0, p2 = 1.0, sr1 = 0.4),
                  resid_enrich = c(all = 0.5),
                  qhat_opt = c(m1 = 8, m2 = 13.33, m3 = 5.75, p1 = 2.77, p2 = 0.72, sr1 = 8.3),    
                  T_opt = c(m1 = 313, m2 = 313, m3 = 303, p1 = 293, p2 = 283, sr1 = 313),
                  T_min = c(m1 = 295.31, m2 = 295.31, m3 = 285.31, p1 = 275.31, p2 = 265.31, sr1 = 273),
                  T_max = c(m1 = 320.67, m2 = 320.67, m3 = 310.67, p1 = 300.67, p2 = 290.67, sr1 = 320.67),
                  ki_NH3_min = c(m1 = 0.01, m2 = 0.015, m3 = 0.015, p1 = 0.015, p2 = 0.015, sr1 = 0.015),
                  ki_NH3_max = c(m1 = 0.10, m2 = 0.131, m3 = 0.131, p1 = 0.131, p2 = 0.131, sr1 = 0.131),
                  ki_NH4_min = c(m1 = 1.70, m2 = 2.714, m3 = 2.714, p1 = 2.714, p2 = 2.714, sr1 = 2.714),
                  ki_NH4_max = c(m1 = 3.10, m2 = 4.764, m3 = 4.764, p1 = 4.764, p2 = 4.764, sr1 = 4.764),
                  pH_upr = c(m1 = 8.0, m2 = 8.0, m3 = 8.0, p1 = 8.0, p2 = 8.0, sr1 = 8.0),
                  pH_lwr = c(m1 = 6.5, m2 = 6.0, m3 = 6.5, p1 = 6.5, p2 = 6.5, sr1 = 6.0)),
  mic_pars = list(ks_SO4 = 0.0067,
                  ki_H2S_meth = 0.23,
                  ki_H2S_sr = 0.25,
                  alpha_opt = 0.015,
                  alpha_T_opt = 313,
                  alpha_T_min = 273,
                  alpha_T_max = 320.67),
  chem_pars = list(COD_conv = c(CH4 = 0.2507, S = 0.5015, VS = 0.69, CO2_anaer = 0.57, CO2_aer = 1.3, CO2_sr = 1.3), kl = c(H2S = 0.032, oxygen = 0.415)),  
  add_pars = NULL,
  startup = -Inf,
  starting = NULL,
  approx_method_temp = 'linear',
  approx_method_pH = 'linear',
  approx_method_SO4 = 'linear',
  par_key = '\\.',
  value = 'ts',                              # Type of output
  warn = TRUE
  ) {

  # If starting conditions are provided from a previous simulation, move to pars
  if (!is.null(starting) & is.data.frame(starting)) {
    message('Using starting conditions from `starting` argument')
    grp_pars[['xa_init']] <- starting[nrow(starting), paste0(names(grp_pars[['qhat_opt']]), '_conc')]
    names(grp_pars[['xa_init']]) <- names(grp_pars[['qhat_opt']])
  }

  # Combine pars to make extraction and pass to rates() easier
  pars <- c(mng_pars, man_pars, grp_pars, mic_pars, chem_pars)

  # Sort out parameter inputs
  # Note: pe.pars = add_pars that use par.element approach, these are converted to normal (simple) add_par elements here
  # Note: sa.pars = normal (simple) add_pars that do not need to be converted
  # Note: Use of [] vs. [[]] affect how code works--needs to work for both lists and vector pars
  # Note: par.element approach is only designed to work for vector elements
  if (!is.null(add_pars) && length(add_pars) > 0 && any(ii <- grepl(par_key, names(add_pars)))) {
    pe.pars <- add_pars[ii]
    sa.pars <- add_pars[!ii]
    apnames <- names(pe.pars)
    pe.pars <- as.numeric(pe.pars)
    split.pars <- strsplit(apnames, par_key)
    pnames <- sapply(split.pars, '[[', 1)
    enames <- sapply(split.pars, '[[', 2)
    names(pe.pars) <- enames
    pe.pars <- split(pe.pars, pnames)
    add_pars <- c(sa.pars, pe.pars)
  }

  # If any additional parameters were added (or modified) using add_pars, update them in pars list here
  if (!is.null(add_pars) && length(add_pars) > 0) {
    if (any(bad.names <- !names(add_pars) %in% names(pars))) {
      stop ('Some `add_pars` names not recognized as valid parameters: ', names(add_pars)[bad.names]) 
    }
    for (i in names(add_pars)) {
      if (!is.data.frame(add_pars[[i]]) && length(pars[[i]]) > 1) {
        pars[[i]][names(add_pars[[i]])] <- add_pars[[i]]
      } else {
        pars[[i]] <- add_pars[[i]]
      }
    }
  }


  # Fill in default values for grp_pars if keyword name `default` or `all` is used
  # Note: Microbial groups are defined by elements in qhat_opt
  # Note: `default1 does *not* work with add_pars argument
  grp_nms <- names(pars$qhat_opt)
  for (i in names(grp_pars)) {
    ppo <- pars[[i]]
    p_nms <- names(pars[[i]])
    if (any(p_nms == 'default')) {
      pars[[i]][grp_nms] <- pars[[i]]['default']
      if (any(p_nms != 'default')) {
        pars[[i]][p_nms[p_nms != 'default']] <- ppo[p_nms[p_nms != 'default']]
      }
    }
    if (any(p_nms == 'all')) {
      pars[[i]][grp_nms] <- pars[[i]]['all']
    }
     # Fix order, drop default element if present
    pars[[i]] <- pars[[i]][grp_nms]
  }

  # Check arguments, including order of element names in some pars
  if (!all.equal(names(pars$yield), names(pars$xa_fresh), names(pars$xa_init), names(pars$decay_rate),
                 names(pars$ks_coefficient), names(pars$resid_enrich), names(pars$qhat_opt), 
                 names(pars$T_opt), names(pars$T_min), names(pars$T_max), 
                 names(pars$ki_NH3_min), names(pars$ki_NH3_max), 
                 names(pars$ki_NH4_min), names(pars$ki_NH4_max), 
                 names(pars$pH_lwr), names(pars$upr))) {
    stop('Microbial groups, i.e., element names in `grp_pars`, must match.')
  }

  # Create temperature function f(t) to allow for variable temperature
  if (is.data.frame(pars$temp_C)) {
    temp <- pars$temp_C$temp_C
    ttime <- pars$temp_C$time
    temp_C_fun <<- approxfun(ttime, temp, method = approx_method_temp, 
                            yleft = temp[1], yright = temp[length(temp)], rule = 2)
  } else {
    temp_C_fun <<- function(x) return(pars$temp_C)
  }

  # Create pH function f(t) to allow for variable pH
  if (is.data.frame(pars$pH)) {
    tpH <- pars$pH$pH
    ttime <- pars$pH$time
    pH_fun <<- approxfun(ttime, tpH, method = approx_method_pH, 
                            yleft = tpH[1], yright = tpH[length(tpH)], rule = 2)
  } else {
    pH_fun <<- function(x) return(pars$pH)
  }

  # Create SO4-2 f(t) for variable acidification with H2SO4
  if (is.data.frame(pars$conc_fresh[['SO4']])) {
    tSO4 <- pars$conc_fres$SO4$SO4
    ttime <- pars$conc_fres$SO4$time
    SO4_fun <<- approxfun(ttime, tSO4, method = approx_method_SO4, 
                            yleft = tSO4[1], yright = tSO4[length(tSO4)], rule = 2)
  } else {
    SO4_fun <<- function(x) return(pars$conc_fresh[['SO4']])
  }

  # Figure out type of run - with constant rates or not
  if (is.numeric(pars$slurry_mass)) {
    # Option 1: Fixed slurry production rate, regular emptying schedule
    dat <- atm_regular(days = days, delta_t = delta_t, pars = pars)
  } else if (is.data.frame(pars$slurry_mass)) {
    # Option 2: Everything based on given slurry mass vs. time
    dat <- atm_variable(days = days, delta_t = delta_t, pars = pars, warn = warn)
  } 

  # Caculate concentrations where relevant
  #conc.names <- names(dat)[!grepl('conc|time|slurry_mass|inhib|qhat|CH4_emis_cum', names(dat))]
  mic_names <- names(pars$qhat)
  conc.names <-  c('NH4', 'NH3', 'Sp', 'VFA', 'sulfide', 'sulfate', mic_names)
  dat_conc <- dat[, conc.names]/dat$slurry_mass
  names(dat_conc) <- paste0(names(dat_conc), '_conc')
  dat <- cbind(dat, dat_conc)

  # Add temperature and pH
  dat$temp_C <- temp_C_fun(dat$time)
  if (is.numeric(pars$pH) | is.data.frame(pars$pH)) {
    dat$pH <- pH_fun(dat$time)
  } else if (pars$pH == 'calc'){
    dat$pH <- H2SO4_titrat(dat$sulfate_conc)
  } else {
    stop('Problem with pH input (bee721)')
  }

  # Add fresh concentrations
  conc_fresh <- unlist(pars$conc_fresh)
  names(conc_fresh) <- paste0(names(conc_fresh), '_conc_fresh')
  dat <- cbind(dat, t(conc_fresh))

  # Calculate rates etc. for ouput, from state variables
  dat$rCH4 <- c(0, diff(dat$CH4_emis_cum))/c(1, diff(dat$time))
  dat$CH4_emis_rate <- c(0, diff(dat$CH4_emis_cum))/c(1, diff(dat$time))
  dat$CH4_emis_rate_slurry <- dat$CH4_emis_rate / (dat$slurry_mass / 1000)
  dat$CH4_flux <- dat$CH4_emis_rate / pars$area
  dat$CO2_emis_rate <- c(0, diff(dat$CO2_emis_cum))/c(1, diff(dat$time))
  dat$CO2_emis_rate_slurry <- dat$CO2_emis_rate / (dat$slurry_mass / 1000)
  dat$CO2_flux <- dat$CO2_emis_rate / pars$area
  ## NTS: Add others, e.g., mu

  # Calculate COD/VS flows
  # First concentrations in g/kg
  dat$dCOD_conc_fresh <- pars$conc_fresh$VFA + pars$conc_fresh$Sp + sum(pars$xa_fresh)
  dat$COD_conc_fresh <- pars$conc_fresh$COD
  dat$ndCOD_conc_fresh <- dat$COD_conc_fresh - dat$dCOD_conc_fresh
  # ndCOD is conserved, same everywhere always
  dat$ndCOD_conc <- ndCOD_conc <- pars$conc_fresh$COD - dat$dCOD_conc_fresh 
  dat$dCOD_conc <- dCOD_conc <- dat$Sp_conc + dat$VFA_conc + rowSums(dat[, paste0(mic_names, '_', 'conc')])
  dat$COD_conc <- COD_conc <- ndCOD_conc + dCOD_conc
  dat$VS_conc <- pars$COD_conv[['VS']] * COD_conc
  # And flows in g/d
  dat$COD_load_rate <- dat$COD_conc_fresh * dat$slurry_prod_rate
  dat$dCOD_load_rate <- dat$dCOD_conc_fresh * dat$slurry_prod_rate
  dat$ndCOD_load_rate <- dat$ndCOD_conc_fresh * dat$slurry_prod_rate
  dat$VS_load_rate <- pars$COD_conv[['VS']] * dat$COD_load_rate
  # Cumulative flow in g
  dat$COD_load_cum <- cumsum(dat$COD_load_rate * delta_t)
  dat$dCOD_load_cum <- cumsum(dat$dCOD_load_rate * delta_t)
  dat$ndCOD_load_cum <- cumsum(dat$ndCOD_load_rate * delta_t)
  dat$VS_load_cum <- cumsum(dat$VS_load_rate * delta_t)
  # And relative emission
  # g CH4/g COD in
  dat$CH4_emis_rate_COD <- dat$CH4_emis_rate / dat$COD_load_rate
  dat$CH4_emis_rate_dCOD <- dat$CH4_emis_rate / dat$dCOD_load_rate
  dat$CH4_emis_rate_VS <- dat$CH4_emis_rate / dat$VS_load_rate
  dat$CH4_emis_cum_COD <- dat$CH4_emis_cum / dat$COD_load_cum
  dat$CH4_emis_cum_dCOD <- dat$CH4_emis_cum / dat$dCOD_load_cum
  dat$CH4_emis_cum_VS <- dat$CH4_emis_cum / dat$VS_load_cum
  # Same for CO2
  dat$CO2_emis_rate_COD <- dat$CO2_emis_rate / dat$COD_load_rate
  dat$CO2_emis_rate_dCOD <- dat$CO2_emis_rate / dat$dCOD_load_rate
  dat$CO2_emis_rate_VS <- dat$CO2_emis_rate / dat$VS_load_rate
  dat$CO2_emis_cum_COD <- dat$CO2_emis_cum / dat$COD_load_cum
  dat$CO2_emis_cum_dCOD <- dat$CO2_emis_cum / dat$dCOD_load_cum
  dat$CO2_emis_cum_VS <- dat$CO2_emis_cum / dat$VS_load_cum
  # Apparent COD conversion fraction
  dat$f_COD_CH4_rate <-  dat$CH4_emis_rate / pars$COD_conv[['CH4']] / dat$COD_load_rate
  #dat$f_COD_CH4_cum <-  dat$CH4_emis_cum / pars$COD_conv[['CH4']] / dat$COD_load_cum
  dat$f_COD_CH4_cum <-  dat$COD_conv_cum_meth / dat$COD_load_cum
  dat$f_COD_respir_cum <-  dat$COD_conv_cum_respir / dat$COD_load_cum
  dat$f_COD_sr_cum <-  dat$COD_conv_cum_sr / dat$COD_load_cum

  # Replace . in names with _
  names(dat) <- gsub('\\.', '_', names(dat))

  # Get averages (excluding startup period)
  dat_sel <- dat[dat$time > startup, ]
  CH4_emis_cum <- dat_sel$CH4_emis_cum[nrow(dat_sel)] - dat_sel$CH4_emis_cum[1]
  CH4_emis_rate <- CH4_emis_cum / (dat_sel$time[nrow(dat_sel)] - dat_sel$time[1])
  # WIP NTS
  CH4_emis_cum

  # Return results
  # Average only
  if (value == 'sum') return(rCH4_ave)
  # ts = time series
  if (value == 'ts') return(dat)
  # tsel = time series after startup period
  if (value == 'tsel') return(dat_sel)
  # Or everything
  return(list(pars = pars, ts = dat, tsel = dat_sel, ave = rCH4_ave))

}

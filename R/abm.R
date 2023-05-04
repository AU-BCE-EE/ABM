abm <-
function(
  days = 365,                               # Number of days to run
  delta_t = 1,                              # Time step for output
  mng_pars = list(slurry_prod_rate = 1000,  # kg/d
                 slurry_mass = 0,           # Initial slurry mass (kg)
                 max_slurry_mass = 33333,   # Maximum slurry mass (kg), default makes for 30 d cycle
                 resid_frac = 0.10,         # Residual slurry fraction after emptying
                 area = 11,                 # Based on 3 m depth (m2)
                 temp_C = 20),
  man_pars = list(conc_fresh = list(S2 = 0.0, SO4 = 0.2, TAN = 1.0, VFA = 4.2, Sp = 65, COD = 160),
                 pH = 7), # Note list not c() so SO4 can be data frame
  grp_pars = list(grps = c('m1', 'm2', 'm3', 'm4', 'm5'),
                  yield = c(all = 0.05),
                  xa_fresh = c(default = 0.001, m3 = 0.01),
                  xa_init = c(m1 = 0.01, m2 = 0.005, m3 = 0.005, m4 = 0.005, m5 = 0.001),
                  decay_rate = c(all = 0.02),
                  ks_coefficient = c(all = 1.0),
                  resid_enrich = c(all = 0.0),
                  qhat_opt = c(m1 = 3.6, m2 = 5.6 , m3 = 7.2, m4 = 8, m5 = 8),
                  T_opt = c(m1 = 18, m2 = 28, m3 = 36, m4 = 43.75, m5 = 55),
                  T_min = c(m1 = 0, m2 = 8, m3 = 15, m4 = 26.25, m5 = 30),
                  T_max = c(m1 = 25, m2 = 38, m3 = 45, m4 = 51.25, m5 = 60),
                  ki_NH3_min = c(all = 0.015),
                  ki_NH3_max = c(all = 0.13),
                  ki_NH4_min = c(all = 2.7),
                  ki_NH4_max = c(all = 4.8),
                  pH_upr = c(all = 8.0),
                  pH_lwr = c(all = 6.5)),
  mic_pars = list(ks_SO4 = 0.0067,
                  ki_H2S_meth = 0.23,
                  ki_H2S_sr = 0.25,
                  alpha_opt = 0.02,
                  alpha_T_min = 0,
                  alpha_T_opt = 50,
                  alpha_T_max = 60),
  chem_pars = list(COD_conv = c(CH4 = 0.2507, S = 0.5015, VS = 0.69, CO2_anaer = 0.53, CO2_aer = 1.1, CO2_sr = 1.2), 
                   kl = c(H2S = 0.02, oxygen = 0.5)),  # kl = mass transfer coefficient (liquid phase units) in m/d
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

  # If starting conditions are provided from a previous simulation, move them to pars
  if (!is.null(starting) & is.data.frame(starting)) {
    message('Using starting conditions from `starting` argument')
    grp_pars[['xa_init']] <- as.numeric(starting[nrow(starting), paste0(names(grp_pars[['qhat_opt']]), '_conc')])
    names(grp_pars[['xa_init']]) <- names(grp_pars[['qhat_opt']])

    mng_pars['slurry_mass'] <- starting[nrow(starting), 'slurry_mass']
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
  # Needs to work in a case where default is all but e.g., m1 is given in add_pars (see def stuff below)
  grp_par_nms <- names(grp_pars)
  grp_par_nms <- grp_par_nms[grp_par_nms != 'grps']
  if (!is.null(add_pars) && length(add_pars) > 0) {
    if (any(bad.names <- !names(add_pars) %in% names(pars))) {
      stop ('Some `add_pars` names not recognized as valid parameters: ', names(add_pars)[bad.names]) 
    }
    # Add in pars (or replace existing elements unless it is time series data added)
    for (i in names(add_pars)) {
      if (!is.data.frame(add_pars[[i]]) && length(pars[[i]]) > 1) {
        pars[[i]][names(add_pars[[i]])] <- add_pars[[i]]
      } else {
        def <- pars[[i]]['all']
        pars[[i]] <- add_pars[[i]]
        if (i %in% grp_par_nms) {
          pars[[i]]['default'] <- def
        }
      }
    }
  }

  # Unlike others, grps in add_pars *will* override default vector (i.e., can be used to remove groups)
  if ('grps' %in% names(add_pars)) {
    pars$grps <- add_pars$grps
  }

  # Fill in default values for grp_pars if keyword name `default` or `all` is used
  # Note: Microbial groups are defined by grps element
  # Note: `default` does *not* work with add_pars argument because grps are already defined in defaults
  # Note: But `all` *does* work
  grp_nms <- pars$grps
  for (i in grp_par_nms) {
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
    # Fix order, drop default element if present, drop unused names
    pars[[i]] <- pars[[i]][grp_nms]
    # Check for missing values
    if (any(is.na(pars[[i]]))) stop('Missing grp_pars elements in ', i, '.')
  }

  # Check arguments, including order of element names in some pars
  # After above block, this should be redundant
  if (!all.equal(names(pars$yield), names(pars$xa_fresh), names(pars$xa_init), names(pars$decay_rate),
                 names(pars$ks_coefficient), names(pars$resid_enrich), names(pars$qhat_opt), 
                 names(pars$T_opt), names(pars$T_min), names(pars$T_max), 
                 names(pars$ki_NH3_min), names(pars$ki_NH3_max), 
                 names(pars$ki_NH4_min), names(pars$ki_NH4_max), 
                 names(pars$pH_lwr), names(pars$upr))) {
    stop('Microbial groups, i.e., element names in `grp_pars`, must match.')
  }

  # Convert temperature constants to K if needed
  for (i in which(grepl('T_', names(pars)))) {
    pars[[i]][pars[[i]] < 200] <- pars[[i]][pars[[i]] < 200] + 273.15
  }

  # Create temperature function f(t) to allow for variable temperature
  if (is.data.frame(pars$temp_C)) {
    temp <- pars$temp_C$temp_C
    ttime <- pars$temp_C$time
    temp_C_fun <- approxfun(ttime, temp, method = approx_method_temp, 
                            yleft = temp[1], yright = temp[length(temp)], rule = 2)
  } else {
    temp_C_fun <- function(x) return(pars$temp_C)
  }

  # Create pH function f(t) to allow for variable pH
  if (is.data.frame(pars$pH)) {
    tpH <- pars$pH$pH
    ttime <- pars$pH$time
    pH_fun <- approxfun(ttime, tpH, method = approx_method_pH, 
                            yleft = tpH[1], yright = tpH[length(tpH)], rule = 2)
  } else {
    pH_fun <- function(x) return(pars$pH)
  }

  # Create SO4-2 f(t) for variable acidification with H2SO4
  if (is.data.frame(pars$conc_fresh[['SO4']])) {
    tSO4 <- pars$conc_fres$SO4$SO4
    ttime <- pars$conc_fres$SO4$time
    SO4_fun <- approxfun(ttime, tSO4, method = approx_method_SO4, 
                            yleft = tSO4[1], yright = tSO4[length(tSO4)], rule = 2)
  } else {
    SO4_fun <- function(x) return(pars$conc_fresh[['SO4']])
  }

  # Figure out type of run - with constant rates or not
  if (is.numeric(pars$slurry_mass)) {
    # Option 1: Fixed slurry production rate, regular emptying schedule
    dat <- abm_regular(days = days, delta_t = delta_t, pars = pars, starting = starting, temp_C_fun = temp_C_fun, pH_fun = pH_fun, SO4_fun = SO4_fun)
  } else if (is.data.frame(pars$slurry_mass)) {
    # Option 2: Everything based on given slurry mass vs. time
    dat <- abm_variable(days = days, delta_t = delta_t, pars = pars, warn = warn, temp_C_fun = temp_C_fun, pH_fun = pH_fun, SO4_fun = SO4_fun)
  } 

  # Caculate concentrations where relevant
  #conc.names <- names(dat)[!grepl('conc|time|slurry_mass|inhib|qhat|CH4_emis_cum', names(dat))]
  mic_names <- pars$grps
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

  # Cut startup period before calculating cumulative flows
  dat <- dat[dat$time > startup, ]

  # Subtract cumulative values that already exist (CH4 emission)
  dat[, grepl('_cum', names(dat))] <- dat[, grepl('_cum', names(dat))] - dat[rep(1, nrow(dat)), grepl('_cum', names(dat))]

  # Calculate COD/VS flows
  # First concentrations in g/kg
  dat$dCOD_conc_fresh <- pars$conc_fresh$VFA + pars$conc_fresh$Sp + sum(pars$xa_fresh)
  dat$COD_conc_fresh <- pars$conc_fresh$COD
  dat$ndCOD_conc_fresh <- dat$COD_conc_fresh - dat$dCOD_conc_fresh
  # ndCOD is conserved, same everywhere always
  dat$ndCOD_conc <- ndCOD_conc <- pars$conc_fresh$COD - dat$dCOD_conc_fresh 
  dat$dCOD_conc <- dCOD_conc <- dat$Sp_conc + dat$VFA_conc + rowSums(dat[, paste0(mic_names, '_', 'conc'), drop = FALSE])
  dat$COD_conc <- COD_conc <- ndCOD_conc + dCOD_conc
  dat$VS_conc <- pars$COD_conv[['VS']] * COD_conc
  # And flows in g/d
  dat$COD_load_rate <- dat$COD_conc_fresh * dat$slurry_prod_rate
  dat$dCOD_load_rate <- dat$dCOD_conc_fresh * dat$slurry_prod_rate
  dat$ndCOD_load_rate <- dat$ndCOD_conc_fresh * dat$slurry_prod_rate
  dat$VS_load_rate <- pars$COD_conv[['VS']] * dat$COD_load_rate
  # Cumulative flow in g
  dat$COD_load_cum <- cumsum(dat$COD_load_rate * c(0, diff(dat$time))) + dat$COD_conc[1] * dat$slurry_mass[1]
  dat$dCOD_load_cum <- cumsum(dat$dCOD_load_rate * c(0, diff(dat$time))) + dat$dCOD_conc[1]* dat$slurry_mass[1]
  dat$ndCOD_load_cum <- cumsum(dat$ndCOD_load_rate * c(0, diff(dat$time))) + dat$ndCOD_conc[1]* dat$slurry_mass[1]
  dat$VS_load_cum <- cumsum(dat$VS_load_rate * c(0, diff(dat$time))) + dat$VS_conc[1]* dat$slurry_mass[1]
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
  dat$f_COD_CH4_cum <-  dat$COD_conv_cum_meth / dat$COD_load_cum
  dat$f_COD_respir_cum <-  dat$COD_conv_cum_respir / dat$COD_load_cum
  dat$f_COD_sr_cum <-  dat$COD_conv_cum_sr / dat$COD_load_cum

  # Replace . in names with _
  names(dat) <- gsub('\\.', '_', names(dat))

  COD_load <- dat$COD_load_cum[nrow(dat)]
  dCOD_load <- dat$dCOD_load_cum[nrow(dat)]
  ndCOD_load <- dat$ndCOD_load_cum[nrow(dat)]
  VS_load <- dat$VS_load_cum[nrow(dat)]

  CH4_emis_cum <- dat$CH4_emis_cum[nrow(dat)]
  CH4_emis_rate <- CH4_emis_cum / (dat$time[nrow(dat)] - dat$time[1])
  CH4_emis_COD <- CH4_emis_cum / COD_load
  CH4_emis_dCOD <- CH4_emis_cum / dCOD_load
  CH4_emis_VS <- CH4_emis_cum / VS_load

  CO2_emis_cum <- dat$CO2_emis_cum[nrow(dat)] - dat$CO2_emis_cum[1]
  CO2_emis_rate <- CO2_emis_cum / (dat$time[nrow(dat)] - dat$time[1])
  CO2_emis_COD <- CO2_emis_cum / COD_load
  CO2_emis_dCOD <- CO2_emis_cum / dCOD_load
  CO2_emis_VS <- CO2_emis_cum / VS_load

  COD_conv_meth <- dat$COD_conv_cum_meth[nrow(dat)] - dat$COD_conv_cum_meth[1]
  COD_conv_respir <- dat$COD_conv_cum_respir[nrow(dat)] - dat$COD_conv_cum_respir[1]
  COD_conv_sr <- dat$COD_conv_cum_sr[nrow(dat)] - dat$COD_conv_cum_sr[1]
  f_COD_CH4 <- COD_conv_meth / COD_load
  f_COD_respir <- COD_conv_respir / COD_load
  f_COD_sr <- COD_conv_sr / COD_load

  summ <- c(COD_load = COD_load, dCOD_load = dCOD_load, ndCOD_load = ndCOD_load, VS_load = VS_load,
            CH4_emis_cum = CH4_emis_cum, CH4_emis_rate = CH4_emis_rate, CH4_emis_COD = CH4_emis_COD, CH4_emis_dCOD = CH4_emis_dCOD, CH4_emis_VS = CH4_emis_VS,
            CO2_emis_cum = CO2_emis_cum, CO2_emis_rate = CO2_emis_rate, CO2_emis_COD = CO2_emis_COD, CO2_emis_dCOD = CO2_emis_dCOD, CO2_emis_VS = CO2_emis_VS,
            COD_conv_meth = COD_conv_meth, COD_conv_respir = COD_conv_respir, COD_conv_sr = COD_conv_sr,
            f_COD_CH4 = f_COD_CH4, f_COD_respir = f_COD_respir, f_COD_sr = f_COD_sr) 


  # Return results
  # Average only
  if (substring(value, 1, 3) == 'sum') return(summ)
  # ts = time series
  if (value == 'ts') return(dat)
  # tsel = time series after startup period
  if (value == 'tsel') return(dat_sel)
  # Or everything
  return(list(pars = pars, ts = dat, tsel = dat_sel, summ = summ))

}

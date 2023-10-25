#abm
abm <- function(
  days = 365,                                # Number of days to run
  delta_t = 1,                               # Time step for output
  times = NA,
  wthr_pars = list(temp_air_C = 20, RH = 90, rain = 1.9, pres_kpa = 101, rs = 10),
  evap_pars = list(evap = 0.5 * et(temp_C = wthr_pars$temp_air_C, pres_kpa = wthr_pars$pres_kpa, rs = wthr_pars$rs)),                # mm/d
  mng_pars = list(slurry_prod_rate = 5500,   # kg/d
                  slurry_mass = 39000,           # Initial slurry mass (kg) NTS: convert to depth??
                  storage_depth = 0.6,         # Storge structure depth, assued to be maximum slurry depth (m)
                  resid_depth = 0.05,         # Residual slurry depth after emptying (m)
                  floor_area = 650,           # Currently for NH3 loss from barn floor (nothing to do with pit/tank floor)
                  area = 650,                 # Area (assume vertical sides) (m2)
                  empty_int = 28,            # (d)
                  temp_C = 20,
                  wash_water = 0,            
                  wash_int = NA,
                  rest_d = 0,
                  cover = NA,
                  resid_enrich = 0.9,
                  slopes = c(urea = NA, slurry_prod_rate = NA),
                  scale = c(ks_coefficient = 1, qhat_opt = 1, xa_fresh = 1, yield = 1, alpha_opt = 1)),
  grz_pars = list(graze_start = "may",
                  graze_days = 0,
                  graze_hours = 0),
  man_pars = list(conc_fresh = list(sulfide = 0.01, urea = 2.4, sulfate = 0.2, TAN = 0.6, starch = 0, 
                                    VFA = 2.83, xa_dead = 0, Cfat = 0, CP = 0, RFd = 0, iNDF = 15, VSd = 75, VSd_A = 44.4, VSnd_A = 20, ash = 15), pH = 7, dens = 1000),
  init_pars = list(conc_init = man_pars$conc_fresh),
  grp_pars = list(grps = c('m0', 'm1', 'm2','sr1'),
                  yield = c(default = 0.05, sr1 = 0.065),
                  xa_fresh = c(default = 0.0628),
                  xa_init = c(all = 0.0628),
                  decay_rate = c(all = 0.02),
                  ks_coefficient = c(default = 1.153337, sr1 = 0.461335),
                  qhat_opt = c(m0 = 0.5515672, m1 = 0.5883362, m2 = 1.5692192, m3 = 7.2, m4 = 8, m5 = 8, sr1 = 8),
                  T_opt = c(m0 = 18, m1 = 18, m2 = 28, m3 = 36, m4 = 43.75, m5 = 55, sr1 = 43.75),
                  T_min = c(m0 = 0, m1 = 8.1928624, m2 = 8.1928624, m3 = 15, m4 = 26.25, m5 = 30, sr1 = 0),
                  T_max = c(m0 = 25, m1 = 25, m2 = 38, m3 = 45, m4 = 51.25, m5 = 60, sr1 = 51.25),
                  ki_NH3_min = c(all = 0.015),
                  ki_NH3_max = c(all = 0.13),
                  ki_NH4_min = c(all = 2.7),
                  ki_NH4_max = c(all = 4.8),
                  ki_H2S_slope = c(default = -0.10623, sr1 = -0.1495),
                  ki_H2S_int = c(default = 0.93066, sr1 = 1.2938),
                  ki_H2S_min = c(default = 0.08)),
  mic_pars = list(ks_SO4 = 0.00694,
                  km_urea = 0.913,
                  alpha_opt = c(urea = 60, VSd = 0.04954023),
                  alpha_T_min = c(urea = 0, VSd = 0),
                  alpha_T_opt = c(urea = 50, VSd = 50),
                  alpha_T_max = c(urea = 60, VSd = 60)),
  chem_pars = list(COD_conv = c(CH4 = 1/0.2507, xa_dead = 1/0.73, RFd = 1/0.8444792, iNDF = 1/0.8444792, starch = 1/0.8444792, 
                                Cfat = 1/0.3117844, CP = 1/0.6541602, VFA = 1/0.9383125, S = 1/0.5015, VS = 1/0.69, CO2_anaer = 1/0.53, CO2_aer = 1/0.436, CO2_sr = 1/1.2, CO2_ureo = 1/1.57,
                                N_CP = 1/0.1014, C_xa_dead = 1/0.358, C_RFd = 1/0.376, C_iNDF = 1/0.358,
                                C_starch = 1/0.377, C_Cfat = 1/0.265, C_CP = 1/0.359 , C_VFA = 1/0.374, C_VSd = 1/0.344, C_N_urea = 1/0.429), 
                   kl = c(NH3 = 54, NH3_floor = 23, H2S = 0.02)), 
  arrh_pars = list(lnA = c(VSd_A = 31.3),
                   E_CH4 = c(VSd_A = 81000), 
                   A = c(xa_dead= 8.56*10^7, starch = 5.86*10^18, Cfat = 0, CP = 181.8, RFd = 1.499476 * 10^12), 
                   E = c(xa_dead= 60600, starch = 109400, Cfat = 0, CP = 23890, RFd = 81052),  
                   R = 8.314,  
                   VS_CH4 = 6.67),
  
  # kl = mass transfer coefficient (liquid phase units) in m/d
  add_pars = NULL,
  pars = NULL,
  startup = 0, # Now number of times complete simulation should be run before returning results
  starting = NULL,
  approx_method_temp = 'linear', # NTS: combine these in 1 vector?
  approx_method_pH = 'linear',
  approx_method_vent_air = 'linear',
  par_key = '\\.',
  value = 'ts',                              # Type of output
  warn = TRUE) {
  
  # If startup repetitions are requested, repeat some number of times before returning results
  if (startup > 0) {

    cat('\nRepeating ')

    for (i in 1:(startup + 1)) {
      cat(i, 'x ')

      # Call abm() with arguments given in outside call except for startup and value
      out <- abm(days = days, delta_t = delta_t, wthr_pars = wthr_pars, evap_pars = evap_pars,
                 mng_pars = mng_pars, grz_pars = grz_pars, man_pars = man_pars, init_pars = init_pars,
                 grp_pars = grp_pars, mic_pars = mic_pars, chem_pars = chem_pars, arrh_pars = arrh_pars,
                 add_pars = add_pars, pars = pars, 
                 startup = 0, 
                 starting = starting, 
                 approx_method_temp = approx_method_temp, approx_method_pH = approx_method_pH, 
                 approx_method_vent_air = approx_method_vent_air, 
                 par_key = par_key, value = 'ts', warn = warn)
 
      # Pull starting *concentrations* (inlcuding xa) from previous sim
      tso <- out

      # Names need to deal with possible data frame for conc_fresh
      cf_names <- names(man_pars$conc_fresh)
      cf_names <- cf_names[!grepl('^time', paste0(cf_names, '_conc'))]

      init_pars$conc_init <- unlist(tso[nrow(tso), paste0(cf_names, '_conc')])
      names(init_pars$conc_init) <- cf_names

      grp_pars$xa_init <- unlist(tso[nrow(tso), paste0(grp_pars$grps, '_conc')])
      names(grp_pars$xa_init) <- grp_pars$grps
 
    }

    # Return only final values
    return(out)
  } 

  # If starting conditions are provided from a previous simulation, move them to pars
  # Note that additional state variables are extracted from starting in abm_*.R
  if (!is.null(starting) & is.data.frame(starting)) {
    message('Using starting conditions from `starting` argument')
    grp_pars[['xa_init']] <- as.numeric(starting[nrow(starting), paste0(names(grp_pars[['qhat_opt']]), '_conc')])
    names(grp_pars[['xa_init']]) <- names(grp_pars[['qhat_opt']])
    mng_pars['slurry_mass'] <- starting[nrow(starting), 'slurry_mass']
  }
  
  # if variable conc fresh, we need to modify the conc_init a little
  if (is.data.frame(man_pars$conc_fresh) & (length(init_pars$conc_init) == length(man_pars$conc_fresh))) {
    init_pars <- list(conc_init = man_pars$conc_fresh[1, -which(names(man_pars$conc_fresh) == "time")])
  } 
  
  
  # Combine pars to make extraction and pass to rates() easier
  if (is.null(pars)) {
    pars <- c(wthr_pars, evap_pars, mng_pars, man_pars, init_pars, grp_pars, mic_pars, chem_pars, grz_pars, arrh_pars, list(days = days))
  }

  # Create error if batch time is not determined and slurry rate should increase over a batch
  if (is.na(pars$wash_int) && !is.na(pars$slopes['slurry_prod_rate'])){
    stop('wash interval cannot be "NA" when slurry production rate increase over a batch of animals, since wash interval defines a batch time')
  }
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
  
  # NTS: Add comment on this (FD?)
  
  # Below code is used only when we have variable concentration of methanogens in the fresh slurry.
  # In that case xa_fresh['time'] will need to be defined in a data.frame, but the block starting below this one 
  # will remove xa_fresh["time"], so we save it as xa_fresh_time before this happens. xa_fresh_time is used in L189 
  
  if(is.data.frame(pars$xa_fresh)) xa_fresh_time <- pars$xa_fresh["time"]
  
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
  
  # Check grp arguments, including order of element names in some pars
  # After above block, this should be redundant
  checkGrpNames(pars)
  
  # Convert temperature constants to K if needed
  pars <- tempsC2K(pars, ll = 200)
  
  # Create time-variable functions
  # Note that pars$x must be numeric constant or df with time (col 1) and var (2)
  temp_C_fun <- makeTimeFunc(pars$temp_C, approx_method = approx_method_temp)
  pH_fun <- makeTimeFunc(pars$pH, approx_method = approx_method_pH)
  conc_fresh_fun <- makeConcFunc(pars$conc_fresh)
  
  # add time to xa_fresh and get xa_fresh funs
  if(is.data.frame(pars$xa_fresh)) pars$xa_fresh$time <- xa_fresh_time
  xa_fresh_fun <- makeXaFreshFunc(pars$xa_fresh)
  
  # Inhibition function
  td <- data.frame(SO4_conc = c(0, 0.09, 0.392251223, 0.686439641, 1.046003263, 1.470942088, 1.961256117, 4), 
                   SO4_inhib = c(1, 1, 0.85, 0.3172, 0.29, 0.1192, 0.05, 0.001))
  SO4_inhibition_fun <- approxfun(td$SO4_conc, td$SO4_inhib, rule = 2)
  
  
  # when variable xa_fresh is used xa_init should be the same as starting xa_fresh
  #if (is.data.frame(pars$xa_fresh)) {
  #  pars$xa_init <- as.numeric(pars$xa_fresh[1, grepl("m[0-9]|sr[0-9]", names(pars$xa_fresh))])
  #  names(pars$xa_init) <- names(pars$xa_fresh)[!grepl("time", names(pars$xa_fresh))]
  #  pars$xa_fresh["time"] <- xa_fresh_time
  #}

  # Convert some supplied parameters
  # Maximum slurry mass in kg
  pars$max_slurry_mass <- pars$storage_depth * pars$area * pars$dens
  pars$resid_mass <- pars$resid_depth / pars$storage_depth * pars$max_slurry_mass

  # Cover effect on NH3 emission rate and N2O
  # reduction from cover 
  pars$EF_NH3 <- coverfun(pars$cover)
  pars$EF_N2O <- ifelse(is.na(pars$cover) | pars$cover == 'none', 0, ifelse(pars$cover == 'tent', 0.05093388, 0.2546694)) # from D. S. Chianese, C. A. Rotz, T. L. Richard, 2009
  
  # calculate grazing interval of year if needed
  if(pars$graze_days > 0){
    pars$graze_int <- c(doy(pars$graze_start, 'March')$day, doy(pars$graze_start, 'March')$day + pars$graze_days)
  } else {
    pars$graze_int <- 0
  }

  if(is.data.frame(pars$slurry_mass)){
    # If missing, set initial slurry mass to 0
    if (pars$slurry_mass[1, 'time'] > 0) {
      pars$slurry_mass <- rbind(c(0, 0), pars$slurry_mass)
    }
    slurry_mass_init <- pars$slurry_mass[1, 'slurry_mass']
  } else {
    slurry_mass_init <- pars$slurry_mass
  }
  
  if (slurry_mass_init == 0) {
    slurry_mass_init <- 1E-10
  }
  
  y <- c(xa = pars$xa_init * slurry_mass_init,
         slurry_mass = slurry_mass_init, 
         xa_dead = pars$conc_init[['xa_dead']] * slurry_mass_init, 
         RFd = pars$conc_init[['RFd']] * slurry_mass_init,
         iNDF = pars$conc_init[['iNDF']] * slurry_mass_init,
         ash = pars$conc_init[['ash']] * slurry_mass_init,
         VSd = pars$conc_init[['VSd']] * slurry_mass_init,
         starch = pars$conc_init[['starch']] * slurry_mass_init,
         CP = pars$conc_init[['CP']] * slurry_mass_init,
         Cfat = pars$conc_init[['Cfat']] * slurry_mass_init,
         VFA = pars$conc_init[['VFA']] * slurry_mass_init, 
         urea = pars$conc_init[['urea']] * slurry_mass_init, 
         TAN = pars$conc_init[['TAN']] * slurry_mass_init, 
         sulfate = pars$conc_init[['sulfate']] * slurry_mass_init, 
         sulfide = pars$conc_init[['sulfide']] * slurry_mass_init, 
         VSd_A = pars$conc_init[['VSd_A']] * slurry_mass_init,
         VSnd_A = pars$conc_init[['VSnd_A']] * slurry_mass_init,
         CH4_A_emis_cum = 0,
         NH3_emis_cum = 0, N2O_emis_cum = 0, CH4_emis_cum = 0, CO2_emis_cum = 0, 
         COD_conv_cum = 0, COD_conv_cum_meth = 0, 
         COD_conv_cum_respir = 0, COD_conv_cum_sr = 0)
  
  if (!is.null(starting) & is.data.frame(starting)) {
    start.vars <- c('slurry_mass', 'xa_dead', 'iNDF', 'ash', 'RFd', 'VSd', 'starch', 'CP', 'Cfat', 'VFA', 'urea', 'TAN', 'sulfate', 'sulfide', 'VSd_A', 'VSnd_A')
    y[start.vars]  <- starting[nrow(starting), start.vars]
  }  
  
  if (any(pars$conc_fresh[['VSd']] > 1e-10) & any(pars$conc_fresh[names(pars$conc_fresh) %in% names(pars$A)] > 0)) {
    stop('Cannot have both VSd and other organic matter components being above 0')
  }
  
  if (is.numeric(pars$slurry_mass)) {
    # Option 1: Fixed slurry production rate, regular emptying schedule
    dat <- abm_regular(days = days, delta_t = delta_t, y = y, pars = pars, starting = starting, temp_C_fun = temp_C_fun, pH_fun = pH_fun,  
                       SO4_inhibition_fun = SO4_inhibition_fun, conc_fresh_fun = conc_fresh_fun, xa_fresh_fun = xa_fresh_fun)
  } else if (is.data.frame(pars$slurry_mass)) {
    # Option 2: Everything based on given slurry mass vs. time
    dat <- abm_variable(days = days, delta_t = delta_t, times = times, y = y, pars = pars, warn = warn, temp_C_fun = temp_C_fun, pH_fun = pH_fun, 
                        SO4_inhibition_fun = SO4_inhibition_fun, conc_fresh_fun = conc_fresh_fun, xa_fresh_fun = xa_fresh_fun)
  } 
  
  colnames(dat) <- gsub("conc_fresh.","conc_fresh_", colnames(dat))
  
  # Calculate concentrations where relevant
  #conc.names <- names(dat)[!grepl('conc|time|slurry_mass|inhib|qhat|CH4_emis_cum', names(dat))]
  mic_names <- pars$grps
  eff_names <- names(dat[grepl("_eff$", names(dat))])
  eff_conc_names <- eff_names[eff_names != "slurry_mass_eff"]
  conc_names <-  c('TAN', 'xa_dead', 'urea', 'RFd', 'iNDF', 'ash', 'VSd', 'starch', 'Cfat', 'CP', 'VFA', 'sulfide', 'sulfate', 'VSd_A', 'VSnd_A', mic_names)
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
  dat$rNH3 <- c(0, diff(dat$NH3_emis_cum))/c(1, diff(dat$time))
  dat$NH3_emis_rate <- c(0, diff(dat$NH3_emis_cum))/c(1, diff(dat$time))
  dat$NH3_emis_rate_slurry <- dat$NH3_emis_rate / (dat$slurry_mass / 1000)
  dat$NH3_flux <- dat$NH3_emis_rate / pars$area
  dat$rN2O <- c(0, diff(dat$N2O_emis_cum))/c(1, diff(dat$time))
  dat$N2O_emis_rate <- c(0, diff(dat$N2O_emis_cum))/c(1, diff(dat$time))
  dat$N2O_emis_rate_slurry <- dat$N2O_emis_rate / (dat$slurry_mass / 1000)
  dat$N2O_flux <- dat$N2O_emis_rate / pars$area
  dat$rCH4 <- c(0, diff(dat$CH4_emis_cum))/c(1, diff(dat$time))
  dat$rCH4_A <- c(0, diff(dat$CH4_A_emis_cum))/c(1, diff(dat$time))
  dat$CH4_emis_rate <- c(0, diff(dat$CH4_emis_cum))/c(1, diff(dat$time))
  dat$CH4_emis_rate_slurry <- dat$CH4_emis_rate / (dat$slurry_mass / 1000)
  dat$CH4_A_emis_rate <- c(0, diff(dat$CH4_A_emis_cum))/c(1, diff(dat$time))
  dat$CH4_flux <- dat$CH4_emis_rate / pars$area
  dat$CO2_emis_rate <- c(0, diff(dat$CO2_emis_cum))/c(1, diff(dat$time))
  dat$CO2_emis_rate_slurry <- dat$CO2_emis_rate / (dat$slurry_mass / 1000)
  dat$CO2_flux <- dat$CO2_emis_rate / pars$area
  
  ## NTS: Add others, e.g., mu
  
  # Calculate COD/VS flows
  # First concentrations in g/kg
  dat$dCOD_conc_fresh <- dat$conc_fresh_VFA + dat$conc_fresh_xa_dead + dat$conc_fresh_RFd + dat$conc_fresh_starch + dat$conc_fresh_CP + dat$conc_fresh_Cfat + dat$conc_fresh_VSd + rowSums(dat[, grepl("xa_fresh_", colnames(dat))])
  dat$COD_conc_fresh <- dat$conc_fresh_VFA + dat$conc_fresh_xa_dead + dat$conc_fresh_RFd + dat$conc_fresh_starch + dat$conc_fresh_CP + dat$conc_fresh_Cfat + dat$conc_fresh_VSd + rowSums(dat[, grepl("xa_fresh_", colnames(dat))]) + dat$conc_fresh_iNDF
  dat$ndCOD_conc_fresh <- dat$COD_conc_fresh - dat$dCOD_conc_fresh
  dat$C_conc_fresh <- dat$conc_fresh_VFA / pars$COD_conv[['C_VFA']] + dat$conc_fresh_xa_dead / pars$COD_conv[['C_xa_dead']] + dat$conc_fresh_RFd / pars$COD_conv[["C_RFd"]] +
                      dat$conc_fresh_starch / pars$COD_conv[['C_starch']] + dat$conc_fresh_CP / pars$COD_conv[['C_CP']] + dat$conc_fresh_Cfat / pars$COD_conv[['C_Cfat']] +
                      dat$conc_fresh_VSd / pars$COD_conv[['C_VSd']] + dat$conc_fresh_iNDF / pars$COD_conv[['C_iNDF']] + rowSums(dat[, grepl("xa_fresh_", colnames(dat))]) / pars$COD_conv[['C_xa_dead']] +
                      dat$conc_fresh_urea / pars$COD_conv[['C_N_urea']]
  
  # ndCOD is almost conserved, same everywhere always
  dat$ndCOD_conc <- ndCOD_conc <- dat$COD_conc_fresh - dat$dCOD_conc_fresh 
  dat$dCOD_conc <- dCOD_conc <- dat$xa_dead_conc + dat$RFd_conc + dat$Cfat_conc + dat$CP_conc + dat$starch_conc + dat$VFA_conc + dat$VSd_conc + rowSums(dat[, paste0(mic_names, '_', 'conc'), drop = FALSE])
  dat$COD_conc <- COD_conc <- ndCOD_conc + dCOD_conc
  dat$VS_conc <- COD_conc / pars$COD_conv[['VS']] # this needs to be updated. Note this is both deg and undeg VS, whereas if VSd is given it is only deg VS
  dat$C_conc <- dat$VFA_conc / pars$COD_conv[['C_VFA']] + dat$xa_dead_conc / pars$COD_conv[['C_xa_dead']] + dat$RFd_conc / pars$COD_conv[["C_RFd"]] +
                dat$starch_conc / pars$COD_conv[['C_starch']] + dat$CP_conc / pars$COD_conv[['C_CP']] + dat$Cfat_conc / pars$COD_conv[['C_Cfat']] +
                dat$VSd_conc / pars$COD_conv[['C_VSd']] + dat$iNDF_conc / pars$COD_conv[['C_iNDF']] + rowSums(dat[, paste0(mic_names, '_', 'conc'), drop = FALSE]) / pars$COD_conv[['C_xa_dead']] +
                dat$urea_conc / pars$COD_conv[['C_N_urea']]
  dat$C_eff_conc <- dat$VFA_eff_conc / pars$COD_conv[['C_VFA']] + dat$xa_dead_eff_conc / pars$COD_conv[['C_xa_dead']] + dat$RFd_eff_conc / pars$COD_conv[["C_RFd"]] +
                dat$starch_eff_conc / pars$COD_conv[['C_starch']] + dat$CP_eff_conc / pars$COD_conv[['C_CP']] + dat$Cfat_eff_conc / pars$COD_conv[['C_Cfat']] +
                dat$VSd_eff_conc / pars$COD_conv[['C_VSd']] + dat$iNDF_eff_conc / pars$COD_conv[['C_iNDF']] + rowSums(dat[, paste0(mic_names, '_', 'conc'), drop = FALSE]) / pars$COD_conv[['C_xa_dead']] +
                dat$urea_eff_conc / pars$COD_conv[['C_N_urea']]
  
  
  dat$Ninorg_conc_fresh <- dat$conc_fresh_urea + dat$conc_fresh_TAN
  dat$Norg_conc_fresh <- dat$conc_fresh_CP / pars$COD_conv[['N_CP']]
  dat$N_conc_fresh <- dat$Ninorg_conc_fresh + dat$Norg_conc_fresh
  dat$Ninorg_conc <- Ninorg_conc <- dat$urea_conc + dat$TAN_conc
  dat$Norg_conc <- Norg_conc <- dat$CP_conc / pars$COD_conv[['N_CP']]
  dat$N_conc <- N_conc <- dat$Ninorg_conc + dat$Norg_conc
  dat$Ninorg_eff_conc <- Ninorg_eff_conc <- dat$urea_eff_conc + dat$TAN_eff_conc
  dat$N_eff_conc <- Ninorg_eff_conc + dat$CP_eff_conc / pars$COD_conv[['N_CP']]

  
  # And flows in g/d
  dat$COD_load_rate <- dat$COD_conc_fresh * dat$slurry_prod_rate
  dat$dCOD_load_rate <- dat$dCOD_conc_fresh * dat$slurry_prod_rate
  dat$ndCOD_load_rate <- dat$ndCOD_conc_fresh * dat$slurry_prod_rate
  dat$VS_load_rate <- dat$COD_load_rate / pars$COD_conv[['VS']]
  dat$C_load_rate <- dat$C_conc_fresh * dat$slurry_prod_rate
  
  dat$Ninorg_load_rate <- dat$Ninorg_conc_fresh * dat$slurry_prod_rate
  dat$Norg_load_rate <- dat$Norg_conc_fresh * dat$slurry_prod_rate
  dat$N_load_rate <- dat$N_conc_fresh * dat$slurry_prod_rate
  
  # Cumulative flow in g
  dat$COD_load_cum <- cumsum(dat$COD_load_rate * c(0, diff(dat$time))) + dat$COD_conc[1] * dat$slurry_mass[1]
  dat$dCOD_load_cum <- cumsum(dat$dCOD_load_rate * c(0, diff(dat$time))) + dat$dCOD_conc[1]* dat$slurry_mass[1]
  dat$ndCOD_load_cum <- cumsum(dat$ndCOD_load_rate * c(0, diff(dat$time))) + dat$ndCOD_conc[1]* dat$slurry_mass[1]
  dat$VS_load_cum <- cumsum(dat$VS_load_rate * c(0, diff(dat$time))) + dat$VS_conc[1]* dat$slurry_mass[1]
  dat$C_load_cum <- cumsum(dat$C_load_rate * c(0, diff(dat$time))) + dat$C_conc[1] * dat$slurry_mass[1]
  
  dat$Ninorg_load_cum <- cumsum(dat$Ninorg_load_rate * c(0, diff(dat$time))) + dat$Ninorg_conc[1] * dat$slurry_mass[1]
  dat$Norg_load_cum <- cumsum(dat$Norg_load_rate * c(0, diff(dat$time))) + dat$Norg_conc[1] * dat$slurry_mass[1]
  dat$N_load_cum <- cumsum(dat$N_load_rate * c(0, diff(dat$time))) + dat$N_conc[1] * dat$slurry_mass[1]

  # And relative emission
  # g CH4/g COD in
  
  dat$CH4_emis_rate_COD <- dat$CH4_emis_rate / dat$COD_load_rate
  dat$CH4_emis_rate_dCOD <- dat$CH4_emis_rate / dat$dCOD_load_rate
  dat$CH4_emis_rate_VS <- dat$CH4_emis_rate / dat$VS_load_rate
  dat$CH4_emis_rate_C <- dat$CH4_emis_rate / dat$C_load_rate
  dat$CH4_emis_cum_COD <- dat$CH4_emis_cum / dat$COD_load_cum
  dat$CH4_emis_cum_dCOD <- dat$CH4_emis_cum / dat$dCOD_load_cum
  dat$CH4_emis_cum_VS <- dat$CH4_emis_cum / dat$VS_load_cum
  dat$CH4_emis_cum_C <- dat$CH4_emis_cum / dat$C_load_cum
  
  # g NH3 N/g N
  dat$NH3_emis_rate_Ninorg <- dat$NH3_emis_rate / dat$Ninorg_load_rate
  dat$NH3_emis_rate_Norg <- dat$NH3_emis_rate / dat$Norg_load_rate
  dat$NH3_emis_rate_N <- dat$NH3_emis_rate / dat$N_load_rate
  dat$NH3_emis_cum_Ninorg <- dat$NH3_emis_cum / dat$Ninorg_load_cum
  dat$NH3_emis_cum_Norg <- dat$NH3_emis_cum / dat$Norg_load_cum
  dat$NH3_emis_cum_N <- dat$NH3_emis_cum / dat$N_load_cum
  
  # g N2O N/g N
  dat$N2O_emis_rate_Ninorg <- dat$N2O_emis_rate / dat$Ninorg_load_rate
  dat$N2O_emis_rate_Norg <- dat$N2O_emis_rate / dat$Norg_load_rate
  dat$N2O_emis_rate_N <- dat$N2O_emis_rate / dat$N_load_rate
  dat$N2O_emis_cum_Ninorg <- dat$N2O_emis_cum / dat$Ninorg_load_cum
  dat$N2O_emis_cum_Norg <- dat$N2O_emis_cum / dat$Norg_load_cum
  dat$N2O_emis_cum_N <- dat$N2O_emis_cum / dat$N_load_cum
  
  # Same for CO2
  dat$CO2_emis_rate_COD <- dat$CO2_emis_rate / dat$COD_load_rate
  dat$CO2_emis_rate_dCOD <- dat$CO2_emis_rate / dat$dCOD_load_rate
  dat$CO2_emis_rate_VS <- dat$CO2_emis_rate / dat$VS_load_rate
  dat$CO2_emis_rate_C <- dat$CO2_emis_rate / dat$C_load_rate
  dat$CO2_emis_cum_COD <- dat$CO2_emis_cum / dat$COD_load_cum
  dat$CO2_emis_cum_dCOD <- dat$CO2_emis_cum / dat$dCOD_load_cum
  dat$CO2_emis_cum_VS <- dat$CO2_emis_cum / dat$VS_load_cum
  dat$CO2_emis_cum_C <- dat$CO2_emis_cum / dat$C_load_cum
  
  # Apparent COD conversion fraction
  dat$f_COD_CH4_rate <-  dat$CH4_emis_rate * pars$COD_conv[['CH4']] / dat$COD_load_rate
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
  CH4_A_emis_cum <- dat$CH4_A_emis_cum[nrow(dat)]
  CH4_emis_rate <- CH4_emis_cum / (dat$time[nrow(dat)] - dat$time[1])
  CH4_A_emis_rate <- CH4_A_emis_cum / (dat$time[nrow(dat)] - dat$time[1])
  CH4_emis_COD <- CH4_emis_cum / COD_load
  CH4_emis_dCOD <- CH4_emis_cum / dCOD_load
  CH4_emis_VS <- CH4_emis_cum / VS_load
  CH4_emis_C <- CH4_emis_cum / C_load
  
  CO2_emis_cum <- dat$CO2_emis_cum[nrow(dat)] - dat$CO2_emis_cum[1]
  CO2_emis_rate <- CO2_emis_cum / (dat$time[nrow(dat)] - dat$time[1])
  CO2_emis_COD <- CO2_emis_cum / COD_load
  CO2_emis_dCOD <- CO2_emis_cum / dCOD_load
  CO2_emis_VS <- CO2_emis_cum / VS_load
  CO2_emis_C <- CO2_emis_cum / C_load
  
  COD_conv_meth <- dat$COD_conv_cum_meth[nrow(dat)] - dat$COD_conv_cum_meth[1]
  COD_conv_respir <- dat$COD_conv_cum_respir[nrow(dat)] - dat$COD_conv_cum_respir[1]
  COD_conv_sr <- dat$COD_conv_cum_sr[nrow(dat)] - dat$COD_conv_cum_sr[1]
  f_COD_CH4 <- COD_conv_meth / COD_load
  f_COD_respir <- COD_conv_respir / COD_load
  f_COD_sr <- COD_conv_sr / COD_load
  
  summ <- c(C_load = C_load, COD_load = COD_load, dCOD_load = dCOD_load, ndCOD_load = ndCOD_load, VS_load = VS_load, Ninorg_load = Ninorg_load, Norg_load = Norg_load, 
            N_load = N_load, NH3_emis_cum = NH3_emis_cum, NH3_emis_rate = NH3_emis_rate, NH3_emis_Ninorg = NH3_emis_Ninorg, NH3_emis_Norg = NH3_emis_Norg,
            NH3_emis_N = NH3_emis_N, N2O_emis_cum = N2O_emis_cum, N2O_emis_rate = N2O_emis_rate, N2O_emis_Ninorg = N2O_emis_Ninorg, N2O_emis_Norg = N2O_emis_Norg,
            N2O_emis_N = N2O_emis_N, CH4_emis_cum = CH4_emis_cum, CH4_emis_rate = CH4_emis_rate, CH4_A_emis_rate = CH4_A_emis_rate, CH4_A_emis_cum = CH4_A_emis_cum, 
            CH4_emis_COD = CH4_emis_COD, CH4_emis_dCOD = CH4_emis_dCOD, CH4_emis_VS = CH4_emis_VS, CH4_emis_C = CH4_emis_C,
            CO2_emis_cum = CO2_emis_cum, CO2_emis_rate = CO2_emis_rate, CO2_emis_COD = CO2_emis_COD, CO2_emis_dCOD = CO2_emis_dCOD, CO2_emis_VS = CO2_emis_VS, CO2_emis_C = CO2_emis_C,
            COD_conv_meth = COD_conv_meth, COD_conv_respir = COD_conv_respir, COD_conv_sr = COD_conv_sr,
            f_COD_CH4 = f_COD_CH4, f_COD_respir = f_COD_respir, f_COD_sr = f_COD_sr) 

  # Add slurry depth
  dat$slurry_depth <- dat$slurry_mass / pars$area / pars$dens

  # Check slurry mass, warn if needed
  if (any(dat$slurry_mass > pars$max_slurry_mass)) {
    warning('Maximum slurry mass exceeded.\nCheck output.')
  }
  
  # Return results
  # Average only
  if (substring(value, 1, 3) == 'sum') return(summ)
  # ts = time series
  if (value == 'ts') return(dat)
  # tsel = time series after startup period
  if (value == 'tsel') return(dat_sel)

  # NTS: create and return cumulative effluent here?
  # NTS: remove duplicate times from emptying (actually only emptying times should show up in effluent right?)
  # NTS: can include in everything below too
  # Or everything
  return(list(pars = pars, ts = dat, tsel = dat_sel, summ = summ))
  
}


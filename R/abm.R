abm <- function(
  days = 365,                               # Number of days to run
  delta_t = 1,                              # Time step for output
  mng_pars = list(slurry_prod_rate = 1000,  # (kg/d)
                  slurry_rem_rate = 0,      # Constant removal rate (kg/d) 
                  slurry_mass = 0,          # Initial slurry mass (kg)
                  storage_depth = 3,        # Storge structure depth, assued to be maximum slurry depth (m)
                  resid_depth = 0.2,        # Residual slurry depth after emptying (m)
                  floor_area = 11,          # NTS: needs to be defined
                  area = 11,                # Area (assume vertical sides) (m2)
                  empty_int = 35,           # (d)
                  temp = 20,                # (deg. C)
                  wash_water = 0,           # (kg) 
                  wash_int = NA,            # (d)
                  rest_d = 0),              # (d)
  man_pars = list(conc_fresh = c(S2 = 0.0, SO4 = 0.2, TAN = 1.0, 
                                 TS = 80, TSS = 60, 
                                 VS = 50, VSS = 30, dsVS = 20, dVSS = 20),
                  pH = 7, dens = 1000), # SO4 cannot be data frame
  #unit_conv = list(conc = '1000 mg/L / g/kg', depth = '0.3048 ft / m', 
  #                 flow = '2.64E-7 mgd / kg/d', temp = '
  #init_pars = list(conc = man_pars$conc_fresh), 
  grp_pars = list(grps = c('m1','m2','m3', 'sr1'),
                  yield = c(default = 0.05, sr1 = 0.065),
                  xa_fresh = c(sr1 = 0, default = 0.02), # (g/kg)?
                  xa_init = c(sr1 = 0, default = 0.02),  # (g/kg)?
                  decay_rate = c(all = 0.02),            # (1/d)
                  ks_coefficient = c(default = 3, sr1 = 1.2),
                  resid_enrich = c(all = 0.1),
                  qhat_opt = c(m1 = 4, m2 = 6, m3 = 8, sr1 = 8),  # (??)
                  T_opt = c(m1 = 20, m2 = 30, m3 = 40, sr1 = 40), # 
                  T_min = c(m1 = 8, m2 = 10, m3 = 20, sr1 = 0),
                  T_max = c(m1 = 25, m2 = 35, m3 = 45, sr1 = 45),
                  ki_NH3_min = c(all = 0.015),
                  ki_NH3_max = c(all = 0.13),
                  ki_NH4_min = c(all = 2.7),
                  ki_NH4_max = c(all = 4.8),
                  pH_upr = c(all = 8.0),
                  pH_lwr = c(all = 6.0)),
  mic_pars = list(ks_SO4 = 0.0067,
                  ki_H2S_meth = 0.23,
                  ki_H2S_sr = 0.25,
                  alpha_opt = 0.01,
                  alpha_T_min = 0,
                  alpha_T_opt = 50,
                  alpha_T_max = 60),
  chem_pars = list(COD_conv = c(CH4 = 0.2507, S = 0.5015, VS = 0.69, CO2_anaer = 0.53, CO2_aer = 1.1, CO2_sr = 1.2), 
                   kl = c(H2S = 0.02, oxygen = 0.5), 
                   unts = list(conc = 'mg/L', depth = 'ft', flow = 'gpm', temp = 'F', mass = 't', area = 'sf')
                  ),  # kl = mass transfer coefficient (liquid phase units) in m/d
  add_pars = NULL,
  startup = -Inf,
  starting = NULL,
  approx_method_temp = 'linear',
  approx_method_pH = 'linear',
  approx_method_SO4 = 'linear',
  par_key = '\\.',
  value = 'ts',                              # Type of output
  nsims = 1,
  warn = TRUE
  ) {

  if (nsims > 1) {
    out <- abm(days = days, delta_t = delta_t, mng_pars = mng_pars, man_pars = man_pars, grp_pars = grp_pars,
               mic_pars = mic_pars, chem_pars = chem_pars, add_pars = add_pars, startup = -Inf, starting = NULL,
               approx_method_temp = approx_method_temp, approx_method_pH = approx_method_pH, approx_method_SO4 = approx_method_SO4, 
               par_key = par_key, value = value, nsims = 1, warn = warn)

    cat('Repeating ')
    for (i in 2:nsims) {
      cat(i, ' ')
      out <- abm(days = days, delta_t = delta_t, mng_pars = mng_pars, man_pars = man_pars, grp_pars = grp_pars,
                 mic_pars = mic_pars, chem_pars = chem_pars, add_pars = add_pars, startup = startup, starting = out,
                 approx_method_temp = approx_method_temp, approx_method_pH = approx_method_pH, approx_method_SO4 = approx_method_SO4, 
                 par_key = par_key, value = value, nsims = 1, warn = warn)
    }

    # Return only final values
    return(out)
  }

  # NTS: fix!
  # NTS: need starting concentrations!
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

  # Unit conversion for inputs and outputs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # NTS: double-check conversion factors
  # Conc to g/kg 
  cf <- list(conc = 1, depth = 1, flow = 1, mass = 1, area = 1)
  if (pars$unts$conc == 'mg/L') {
    cf$conc <- 0.001 / (pars$dens / 1000) # g/kg per mg/L 
  } else if (pars$unts$conc != 'g/kd') {
    stop('Error in input for unts$conc')
  }

  # Depth to m
  if (pars$unts$depth == 'ft') {
    cf$depth <- 0.3048 # m per ft
  } else if (pars$unts$depth != 'm') {
    stop('Error in input for unts$depth')
  }

  # Flow to kg/d
  if (pars$unts$flow == 'mgd') {
    cf$flow <- 1 * (pars$dens / 1000)
  } else if (pars$unts$flow == 'gpd') {
    cf$flow <- 1 * (pars$dens / 1000)
  } else if (pars$unts$flow == 'gpm') {
    cf$flow <- 5451 * (pars$dens / 1000) # kg/d per gal/min
  } else if (pars$unts$flow != 'kg/d') {
    stop('Error in input for unts$flow')
  }

  if (pars$unts$mass == 't') {
    cf$mass <- 1 / 1.102E-3
  } else if (pars$unts$mass != 'kg') {
    stop('Error in input for unts$mass')
  }

  if (pars$unts$area == 'sf') {
    cf$area <- 1 / 10.764
  } else if (pars$unts$area != 'kg') {
    stop('Error in input for unts$area')
  }

  # Apply conversion factors
  pars$conc_fresh <- pars$conc_fresh * cf$conc

  pars$storage_depth <- pars$storage_depth * cf$depth
  pars$resid_depth <- pars$resid_depth * cf$depth 

  pars$slurry_prod_rate <- pars$slurry_prod_rate * cf$flow
  pars$slurry_rem_rate <- pars$slurry_rem_rate * cf$flow

  pars$slurry_mass <- pars$slurry_mass * cf$mass

  pars$floor_area <- pars$floor_area * cf$area
  pars$area <- pars$area * cf$area

  # Temperature trickier
  if (pars$unts$temp == 'F') {
    if (is.data.frame(pars$temp)) {
      pars$temp$temp <- (pars$temp$temp - 32) * 5 / 9
    } else {
      pars$temp <- (pars$temp - 32) * 5 / 9
    }
  } else if (pars$unts$temp != 'C') {
    stop('Error in input for unts$temp')
  }
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


  # Create temperature function f(t) to allow for variable temperature
  if (is.data.frame(pars$temp)) {
    if (!all(names(pars$temp)[1:2] %in% c('time', 'temp'))) {
      stop('Check names in temp data frame')
    }
    temp <- pars$temp$temp
    ttime <- pars$temp$time
    temp_fun <- approxfun(ttime, temp, method = approx_method_temp, 
                            yleft = temp[1], yright = temp[length(temp)], rule = 2)
  } else {
    temp_fun <- function(x) return(pars$temp)
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
  # NTS: no longer an option now with unit conversion!
  # Need to remove!
  if (is.data.frame(pars$conc_fresh[['SO4']])) {
    tSO4 <- pars$conc_fres$SO4$SO4
    ttime <- pars$conc_fres$SO4$time
    SO4_fun <- approxfun(ttime, tSO4, method = approx_method_SO4, 
                            yleft = tSO4[1], yright = tSO4[length(tSO4)], rule = 2)
  } else {
    SO4_fun <- function(x) return(pars$conc_fresh[['SO4']])
  }


  # Calculate fresh concentrations
  # dsVS -> dsOM
  # dVSS -> dpOM
  # Rename a few 
  pars$conc_fresh[['dpVS']] <- pars$conc_fresh[['dVSS']] 
  pars$conc_fresh[['dsVS']] <- pars$conc_fresh[['dsVS']] 
  pars$conc_fresh[['pVS']] <- pars$conc_fresh[['VSS']] 
  # Calculate remainder by difference
  pars$conc_fresh[['sVS']]  <- pars$conc_fresh[['VS']] - pars$conc_fresh[['pVS']]
  pars$conc_fresh[['isVS']] <- pars$conc_fresh[['sVS']] - pars$conc_fresh[['dsVS']] 
  pars$conc_fresh[['ipVS']] <- pars$conc_fresh[['pVS']] - pars$conc_fresh[['dpVS']] 

  pars$conc_fresh[['iFS']] <- pars$conc_fresh[['TS']] - pars$conc_fresh[['VS']] 

  # Convert to COD
  pars$conc_fresh[['COD']]   <- pars$conc_fresh[['VS']]   / pars$COD_conv[['VS']]
  pars$conc_fresh[['pCOD']]  <- pars$conc_fresh[['pVS']]  / pars$COD_conv[['VS']]
  pars$conc_fresh[['dpCOD']] <- pars$conc_fresh[['dpVS']] / pars$COD_conv[['VS']]
  pars$conc_fresh[['ipCOD']] <- pars$conc_fresh[['ipVS']] / pars$COD_conv[['VS']]
  pars$conc_fresh[['sCOD']]  <- pars$conc_fresh[['sVS']]  / pars$COD_conv[['VS']]
  pars$conc_fresh[['dsCOD']] <- pars$conc_fresh[['dsVS']] / pars$COD_conv[['VS']]
  pars$conc_fresh[['isCOD']] <- pars$conc_fresh[['isVS']] / pars$COD_conv[['VS']]

  # Convert some supplied parameters
  # Maximum slurry mass in kg
  pars$max_slurry_mass <- pars$storage_depth * pars$area * pars$dens
  pars$resid_mass <- pars$resid_depth / pars$storage_depth * pars$max_slurry_mass
  
  if(is.data.frame(pars$slurry_mass)){
    slurry_mass_init <- pars$slurry_mass[1, 'slurry_mass']
  } else {
    slurry_mass_init <- pars$slurry_mass
  }
  
  if (slurry_mass_init == 0) {
    slurry_mass_init <- 1E-10
  }
 
  # Initial state variable vector
  y <- c(xa = pars$xa_init, 
         slurry_mass = 1, 
         unlist(pars$conc_fresh[c('dpCOD', 'ipCOD', 'dsCOD', 'isCOD', 'iFS')]), # NTS: unlist prob not needed now that conc_fresh is vector
         sulfate = SO4_fun(0), 
         sulfide = pars$conc_fresh[['S2']]) 
  # Convert to mass (g??)
  y <- y * slurry_mass_init
  # Add 0 for cumulative emission
  zz <- c('CH4_emis_cum', 'CO2_emis_cum', 'COD_conv_cum', 'COD_conv_cum_meth', 
          'COD_conv_cum_respir', 'COD_conv_cum_sr')
  y[zz] <- 0

  # Figure out type of run - with constant rates or not
  if (is.numeric(pars$slurry_mass)) {
    # Option 1: Fixed slurry production rate, regular emptying schedule
    dat <- abm_regular(days = days, delta_t = delta_t, y = y, pars = pars, temp_fun = temp_fun, pH_fun = pH_fun, SO4_fun = SO4_fun)
  } else if (is.data.frame(pars$slurry_mass)) {
    # Option 2: Everything based on given slurry mass vs. time
    dat <- abm_variable(days = days, delta_t = delta_t, y = y, pars = pars, warn = warn, temp_fun = temp_fun, pH_fun = pH_fun, SO4_fun = SO4_fun)
  } 

  # Calculate totals etc.
  dat$COD <- dat$dpCOD + dat$ipCOD + dat$dsCOD + dat$isCOD
  dat$pCOD <- dat$dpCOD + dat$ipCOD
  dat$sCOD <- dat$dsCOD + dat$isCOD
  dat$VS <- pars$COD_conv[['VS']] * dat$COD
  dat$VSS <- pars$COD_conv[['VS']] * dat$pCOD
  dat$TS <- dat$iFS + dat$VS

  # Caculate concentrations where relevant (NTS: g/kg????)
  mic_names <- pars$grps
  conc.names <-  c('NH4', 'NH3', 
                   'COD', 'dpCOD', 'ipCOD', 'dsCOD', 'isCOD', 'pCOD', 'sCOD',
                   'TS', 'VS', 'VSS',
                   'sulfide', 'sulfate', mic_names)
  dat_conc <- dat[, conc.names] / dat$slurry_mass
  names(dat_conc) <- paste0(names(dat_conc), '_conc')
  dat <- cbind(dat, dat_conc)

  # Add temperature and pH
  dat$temp <- temp_fun(dat$time)
  if (is.numeric(pars$pH) | is.data.frame(pars$pH)) {
    dat$pH <- pH_fun(dat$time)
  } else if (pars$pH == 'calc'){
    dat$pH <- H2SO4_titrat(dat$sulfate_conc)
  } else {
    stop('Problem with pH input (bee721)')
  }

  # Add fresh concentrations
  # NTS: Simpler now that input is vector not list
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

  # Calculate COD/VS flows
  # First concentrations in g/kg
  dat$dCOD_conc_fresh <- pars$conc_fresh['dsCOD'] + pars$conc_fresh['dpCOD'] + sum(pars$xa_fresh)
  dat$COD_conc_fresh <- pars$conc_fresh['COD']
  dat$ndCOD_conc_fresh <- dat$COD_conc_fresh - dat$dCOD_conc_fresh

  # Loading 
  # Flows in g/d
  dat$COD_load_rate <- dat$COD_conc_fresh * dat$slurry_prod_rate
  dat$dCOD_load_rate <- dat$dCOD_conc_fresh * dat$slurry_prod_rate
  dat$ndCOD_load_rate <- dat$ndCOD_conc_fresh * dat$slurry_prod_rate
  dat$VS_load_rate <- dat$VS_conc_fresh * dat$slurry_prod_rate

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
  dat$CH4_rate_f_COD <-  dat$CH4_emis_rate / pars$COD_conv[['CH4']] / dat$COD_load_rate
  dat$CH4_cum_f_COD <-  dat$COD_conv_cum_meth / dat$COD_load_cum
  dat$respir_cum_f_COD <-  dat$COD_conv_cum_respir / dat$COD_load_cum
  dat$sr_cum_f_COD <-  dat$COD_conv_cum_sr / dat$COD_load_cum
  dat$loss_cum_f_COD <-  dat$CH4_cum_f_COD + dat$respir_cum_f_COD + dat$sr_cum_f_COD

  # Add slurry depth
  dat$slurry_depth <- dat$slurry_mass / pars$dens / pars$area

  # Convert units back to inputs
  dat[, grep('conc', names(dat))] <- dat[, grep('conc', names(dat))] / cf$conc
  dat[, grep('area', names(dat))] <- dat[, grep('area', names(dat))] / cf$area
  dat[, grep('depth', names(dat))] <- dat[, grep('depth', names(dat))] / cf$depth
  dat[, 'slurry_prod_rate'] <- dat[, 'slurry_prod_rate'] / cf$flow
  dat[, 'slurry_mass'] <- dat[, 'slurry_mass'] / cf$mass

  # Calculate totals for summary
  which.tot <- grep('cum', names(dat), value = TRUE)
  summ <- unlist(dat[nrow(dat), which.tot])

  # Replace . in names with _
  names(dat) <- gsub('\\.', '_', names(dat))

  # Return results
  # Average only
  if (substring(value, 1, 3) == 'sum') return(summ)
  # ts = time series
  if (value == 'ts') return(dat)
  # Or everything
  return(list(pars = pars, ts = dat, summ = summ))

}

# Sorts out parameters and packages them all together in the output list

packPars <- function(mng_pars,
                     man_pars,
                     init_pars,
                     grp_pars,
                     mic_pars,
                     chem_pars,
                     ctrl_pars,
                     add_pars,
                     pars,
                     starting,
                     days) {

  # If starting conditions are provided from a previous simulation, move them to pars
  # Note that additional state variables are extracted from `starting` in abm_*.R
  if (!is.null(starting) & is.data.frame(starting)) {
    message('Using starting conditions from `starting` argument')
    grp_pars[['xa_init']] <- as.numeric(starting[nrow(starting), paste0(names(grp_pars[['qhat_opt']]), '_conc')])
    names(grp_pars[['xa_init']]) <- names(grp_pars[['qhat_opt']])
    mng_pars['slurry_mass'] <- starting[nrow(starting), 'slurry_mass']
  }
  
  # Combine pars to make extraction and pass to rates() easier
  if (is.null(pars)) { 
    #pars <- c(mng_pars, man_pars, init_pars, grp_pars, mic_pars, chem_pars, list(days = days), ctrl_pars)
    pars <- c(mng_pars, man_pars, init_pars, grp_pars, mic_pars, chem_pars, ctrl_pars)
  }

  ## if variable conc fresh, we need to modify the conc_init a little
  #if (is.data.frame(pars$conc_fresh) & (length(pars$conc_init) == length(pars$conc_fresh))) {
  #  pars$conc_init <- pars$conc_fresh[1, -which(names(pars$conc_fresh) == "time")]
  #} 
  
  # Combine pars to make extraction and pass to rates() easier
  # Sort out parameter inputs
  # Note: pe.pars = add_pars that use par.element approach, these are converted to normal (simple) add_par elements here
  # Note: sa.pars = normal (simple) add_pars that do not need to be converted
  # Note: Use of [] vs. [[]] affect how code works--needs to work for both lists and vector pars
  # Note: par.element approach is only designed to work for vector elements
  par_key <- ctrl_pars$par_key
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
  
  ## NTS: Add comment on this (FD?)
  ## Below code is used only when we have variable concentration of methanogens in the fresh slurry.
  ## In that case xa_fresh['time'] will need to be defined in a data.frame, but the block starting below this one 
  ## will remove xa_fresh["time"], so we save it as xa_fresh_time before this happens. xa_fresh_time is used in L189 
  #if(is.data.frame(pars$xa_fresh)) xa_fresh_time <- pars$xa_fresh["time"]
  
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

  # Add number of groups and indices of different groups
  pars$n_mic <- length(pars$grps)
  pars$i_meth <- grep('^[mp]', pars$grps)
  pars$i_sr <- grep('^sr', pars$grps)
  pars$i_mic <- grep('^sr|^p|^m', pars$grps)
  pars$i_hyd <- grep('^hyd', pars$grps)
  pars$i_aer <- grep('^aer', pars$grps)
  
  # Convert temperature constants to K if needed
  pars <- tempsC2K(pars, ll = 200)
  
  # Convert some supplied parameters
  # Maximum slurry mass in kg
  pars$max_slurry_mass <- pars$storage_depth * pars$area * pars$dens
  pars$resid_mass <- pars$resid_depth / pars$storage_depth * pars$max_slurry_mass

  # Cover effect on NH3 emission rate and N2O
  # Reduction from cover 
  pars$EF_NH3 <- coverfun(pars$cover, pars$scale_EF_NH3)
  pars$EF_N2O <- ifelse(pars$cover == 'none', 0, ifelse(pars$cover == 'tent', 0.05093388, 0.2546694)) # from D. S. Chianese, C. A. Rotz, T. L. Richard, 2009
 
  return(pars)
 
}

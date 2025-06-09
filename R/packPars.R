# Sorts out parameters and packages them all together in the output list

packPars <- function(mng_pars,
                     man_pars,
                     init_pars,
                     grp_pars,
                     mic_pars,
                     sub_pars,
                     chem_pars,
                     inhib_pars,
                     ctrl_pars,
                     var_pars,
                     add_pars,
                     pars,
                     starting) {

  # Move extra var_pars
  if(length(var_pars) > 1) {
    for (i in 2:length(var_pars)) {
      dv <- var_pars[[i]]
      nv <- names(var_pars)[[i]]
      if (nrow(dv) != nrow(var_pars$var)) {
        stop('var_pars has multiple elements with different sizes.')
      }
      if (!all(dv$time == var_pars$var$time)) {
        stop('Multiple elements in var_pars have different time values in first column.')
      }
      ll <- list()
      for (j in 1:nrow(var_pars$var)) {
        # Challenging to get list into each element using indexing
        x <- as.numeric(dv[j, -1, drop = FALSE])
        names(x) <- names(dv)[-1]
        ll[[j]] <- x
      }
      var_pars$var[[nv]] <- ll
    }
    # Take only the var element and var_pars must remain a list to avoid duplicate par elements from combining into par
    var_pars <- var_pars['var']
  }

  # Combine pars to make extraction and pass to rates() easier
  if (is.null(pars)) { 
    pars <- c(mng_pars, man_pars, init_pars, grp_pars, mic_pars, sub_pars, chem_pars, inhib_pars, ctrl_pars, var_pars)
  }

  # Sort out more parameter inputs
  # Note: pe.pars = add_pars that use par.element approach, these are converted to normal (simple) add_par elements here
  # Note: sa.pars = normal (simple) add_pars that do not need to be converted
  # Note: Use of [] vs. [[]] affect how code works--needs to work for both lists and vector pars
  # Note: par.element approach is only designed to work for vector elements
  par_key <- ctrl_pars$par_key
  if (!is.null(add_pars) && length(add_pars) > 0 && any(ii <- grepl(par_key, names(add_pars)))) {
    pe.pars <- add_pars[ii]
    sa.pars <- add_pars[!ii]
    apnames <- names(pe.pars)
    pe.pars[!grepl('[A-Za-z]', pe.pars)] <- as.numeric(pe.pars[!grepl('[A-Za-z]', pe.pars)])
    split.pars <- strsplit(apnames, par_key)
    pnames <- sapply(split.pars, '[[', 1)
    enames <- sapply(split.pars, '[[', 2)
    names(pe.pars) <- enames
    pe.pars <- split(pe.pars, pnames)
    add_pars <- c(sa.pars, pe.pars)
  }


  # If any additional parameters were added (or modified) using add_pars, update them in pars list here
  # But grp_pars and sub_pars work differently than the others because of the all = and default = keywords
  # Needs to work in a case where default is all but e.g., m1 is given in add_pars (see def stuff below)
  grp_par_nms <- names(grp_pars)[names(grp_pars) != 'grps']
  sub_par_nms <- names(sub_pars)[names(sub_pars) != 'subs']
  if (!is.null(add_pars) && length(add_pars) > 0) {
    if (any(bad.names <- !names(add_pars) %in% names(pars))) {
      stop ('Some `add_pars` names not recognized as valid parameters: ', names(add_pars)[bad.names]) 
    }
    # Add in pars (or replace existing elements unless it is time series data added)
    for (i in names(add_pars)) {
      if (!is.data.frame(add_pars[[i]]) && length(pars[[i]]) > 1) {
        pars[[i]][names(add_pars[[i]])] <- unlist(add_pars[[i]])
      } else {
        def <- pars[[i]]['all']
        pars[[i]] <- add_pars[[i]]
        if (i %in% c(grp_par_nms, sub_par_nms)) {
          pars[[i]]['default'] <- def
        }
      }
    }
  }
  
  # Unlike others, grps and subs in add_pars *will* override default vector (i.e., can be used to remove groups)
  if ('grps' %in% names(add_pars)) {
    pars$grps <- add_pars$grps
  }
  if ('subs' %in% names(add_pars)) {
    pars$subs <- add_pars$subs
  }

  # Fill in default values for grp_pars if keyword name `default` or `all` is used
  # Note: Microbial groups are defined by grps element
  # Note: `default` does *not* work with add_pars argument because grps are already defined in defaults
  # Note: But `all` *does* work
  # expandPars() will also sort out element order and drop excluded elements
  pars <- expandPars(pars = pars, elnms = pars$grps, parnms = grp_par_nms)
  pars <- expandPars(pars = pars, elnms = pars$subs, parnms = sub_par_nms)

  # Check grp arguments, including order of element names in some pars
  # After above block, this should be redundant
  checkGrpNames(pars)

  # For size-variable parameters, get number of elements and indices
  # NTS: I expect to change to a different approach, only using block below this one with names
  pars$n_mic <- length(pars$grps)
  pars$i_mic <- grep('^sr|^p|^m', pars$grps)
  pars$i_meth <- grep('^[mp]', pars$grps)
  pars$i_sr <- grep('^sr', pars$grps)
  pars$i_aer <- grep('^aer', pars$grps)
  pars$i_hyd <- grep('^hyd', pars$grps)

  # Get names of variable elements
  # Remember pars$grps/pars$mics and pars$subs already exist (set in pars input)
  pars$meths <- pars$grps[pars$i_meth]
  pars$srs <- pars$grps[pars$i_sr]
  pars$aers <- pars$grps[pars$i_aer]
  pars$hyds <- pars$grps[pars$i_hyd]

  # All solutes
  pars$sols <- c(pars$comps, 'VFA')
  pars$conc_fresh <- c(pars$comp_fresh, pars$VFA_fresh)

  # Master species, fill in masters = masters
  mmspec <- pars$sols 
  names(mmspec) <- mmspec
  pars$mspec <- c(pars$mspec, mmspec)
  pars$mspec <- pars$mspec[!duplicated(names(pars$mspec))]
  
  # Substrates
  pars$n_subs <- length(pars$subs)
  pars$i_subs <- length(pars$grps) + 1:length(pars$subs)
  
  # Other constants
  pars$g_NH4 <- 0.7
  pars$temp_standard <- 298
  pars$temp_zero <- 273
  pars$temp_K <- pars$temp_C + 273.15
  pars$pH_floor <- 7

  # Convert temperature constants to K if needed
  pars <- tempsC2K(pars, cutoff = 200)
  
  # Convert some supplied parameters
  # Maximum slurry mass in kg
  pars$max_slurry_mass <- pars$storage_depth * pars$area * pars$dens
  pars$resid_mass <- pars$resid_depth / pars$storage_depth * pars$max_slurry_mass
  
  # If starting conditions are provided from a previous simulation, move them to pars
  # Note that additional state variables are extracted from `starting` in abm_*.R
  if (!is.null(starting) & is.data.frame(starting)) {
    message('Using starting conditions from `starting` argument')
    pars$xa_init[pars$grps] <- as.numeric(starting[nrow(starting), paste0(pars$grps, '_conc')])
    pars$conc_init['VFA'] <- as.numeric(starting[nrow(starting), 'VFA_conc'])
    pars$sub_init[pars$subs] <- as.numeric(starting[nrow(starting), paste0(pars$subs, '_conc')])
    # Set slurry_mass as well?
    # comp solutes??
  }

  return(pars)
 
}

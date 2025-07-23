  
moveStartingPars <- function(pars, starting) {
  
  message('Using starting conditions from `starting` argument')
  pars$xa_init[pars$grps] <- as.numeric(starting[nrow(starting), paste0(pars$grps, '_conc')])
  pars$conc_init['CH3COOH'] <- as.numeric(starting[nrow(starting), 'CH3COOH_conc'])
  pars$sub_init[pars$subs] <- as.numeric(starting[nrow(starting), paste0(pars$subs, '_conc')])
  # Set slurry_mass as well?
  # NTS: Set comp solutes also?

  return(pars)
  
}


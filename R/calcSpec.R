# Calculate chemical speciation
# Always includes VFA
# May include optional parameter-defined species 

calcSpec <- function(pars, y) {

  # Get concentrations of all solutes
  concs <- y[pars$sols] / y['slurry_mass']

  # Save totals 
  pars$conc_tot <- concs
  
  # Calculate and save species concentrations
  pars$conc_sp <- as.numeric(concs[pars$mspec[pars$specs]]) * 1 / (1 + 10^(-pars$lkan[pars$specs] - pars$pH))
  pars$conc_sp[pars$mspec[pars$spec]] <- pars$conc_tot[pars$mspec[pars$spec]] - pars$conc_sp[pars$spec]
  
  return(pars)
}

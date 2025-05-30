# Calculate inhibition from pars values

calcInhib <- function(pars, y) {

  # Get concentrations
  # NTS: Ultimately this needs to reflect full flexibility of abm()
  # NTS: Use function here
  concs <- c(NH4. = pars$conc_fresh[['TAN']], NH3 = 0, VFA = y[['VFA']] / y[['slurry_mass']], HVFA = 0)

  pars$ired <- inhib(pars$ilwr, pars$iupr, concs)

  return(pars)

}

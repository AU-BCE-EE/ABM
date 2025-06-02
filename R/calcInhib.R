# Calculate inhibition from pars values

calcInhib <- function(pars, y) {

  if (is.null(pars$ilwr) | is.null(pars$iupr)) {
    pars$ired <- rep(1, length(pars$grps))
    names(pars$ired) <- pars$grps
    return(pars)
  }

  # Get concentrations
  # NTS: Ultimately this needs to reflect full flexibility of abm()
  # NTS: Use function here
  concs <- c(NH4. = pars$conc_fresh[['TAN']], NH3 = 0, VFA = y[['VFA']] / y[['slurry_mass']], HVFA = 0)

  pars$ired <- inhib(pars$ilwr, pars$iupr, concs)

  return(pars)

}

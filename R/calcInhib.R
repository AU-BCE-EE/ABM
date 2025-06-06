# Calculate inhibition from pars values

calcInhib <- function(pars, y) {
  
  # Skip all this if there are no inhibition parameters
  if (!is.null(pars$ilwr) & !is.null(pars$iupr)) {
    # Get concentrations of inhibitor species
    concs <- pars$conc_sp[colnames(pars$ilwr)]
    pars$ired <- inhib(pars$ilwr, pars$iupr, concs)
  } else {
    pars$ired <- rep(1, length(pars$grps))
    names(pars$ired) <- pars$grps
  }

  return(pars)

}

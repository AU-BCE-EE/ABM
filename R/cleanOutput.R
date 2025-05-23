# Clean up abm_*() output before returning it

cleanOutput <- function(dat, pars, addcols, addconcs) {

  # Replace . in names with _
  names(dat) <- gsub('\\.', '_', names(dat))
  
  if (addcols) {
    # Add slurry depth
    dat$slurry_depth <- dat$slurry_mass / pars$area / pars$dens
  }

  if (addconcs) {
    dat <- addConcs(dat, pars)
  }
  
  return(dat)

}


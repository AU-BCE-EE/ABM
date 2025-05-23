

cleanOutput <- function(dat, pars, addcols) {

  # Replace . in names with _
  names(dat) <- gsub('\\.', '_', names(dat))
  
  if (addcols) {
    # Add slurry depth
    dat$slurry_depth <- dat$slurry_mass / pars$area / pars$dens
  }

  return(dat)

}


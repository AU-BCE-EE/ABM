# Clean up abm_*() output before returning it

cleanOutput <- function(dat, pars, addcols, addconcs, cumeff) {

  # Replace . in names with _
  names(dat) <- gsub('\\.', '_', names(dat))
  
  if (addcols) {
    # Add slurry depth
    dat$slurry_depth <- dat$slurry_mass / pars$area / pars$dens
  }

  # Make effluent results cumulative
  if (cumeff) {
    dat[, names(dat)[grepl('_eff', names(dat))]] <- lapply(dat[, names(dat)[grepl('_eff', names(dat))]], cumsum)
  }

  # Add concentrations (after cumulative effluent)
  if (addconcs) {
    dat <- addConcs(dat, pars)
  }
  
  return(dat)

}


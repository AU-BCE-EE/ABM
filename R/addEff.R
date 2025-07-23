# Add effluent results

addEff <- function(dat, y.eff) {

  dat <- data.frame(dat)

  # Create new effluent columns only if they do not exist
  # Should never be a case when some but not all exist
  if (!all(names(y.eff) %in% names(dat))) {
    dat[, names(y.eff)] <- 0
  }
  
  # Add in effluent results
  dat[nrow(dat), names(y.eff)] <- y.eff

  return(dat)

}

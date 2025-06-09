# Copy time-variable parameters into their normal par position for an interval

updateVarPars <- function(pars, i) {

  if (any(! names(pars$var) %in% c('time', 'slurry_mass'))) {
    vdat <- pars$var[, !grepl('time|slurry_mass', names(pars$var)), drop = FALSE]

    for (j in 1:ncol(vdat)) {
      pn <- names(vdat)[j]
      newval <- vdat[i, j] 
      if (is.list(newval)) {
        newval <- newval[[1]]
      }
      pars[[pn]] <- newval
    }

    # Convert temperature to K in case any C values were given in var
    pars <- tempsC2K(pars, cutoff = 200)
  }

  return(pars)

}

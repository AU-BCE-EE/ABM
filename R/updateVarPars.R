# Copy time-variable parameters into their normal par position for an interval

updateVarPars <- function(pars, y, i) {

  vdat <- pars$var[, names(pars$var) != 'time', drop = FALSE]

  for (j in 1:ncol(vdat)) {
    pn <- names(vdat)[j]
    newval <- vdat[i, j] 
    if (is.list(newval)) {
      newval <- newval[[1]]
    }
    pars[[pn]] <- newval
  }

  # Temperature and temperature-dependent pars could change
  # Convert temperature to K in case any C values were given in var
  pars <- tempsC2K(pars, cutoff = 200)

  # Calculate temperature-dependent par values
  pars <- calcTempPars(pars, y)

  return(pars)

}

# Converts temperature parameters from C to K

tempsC2K <- function(pars, ll) {
  for (i in which(grepl('T_', names(pars)))) {
    pars[[i]][pars[[i]] < ll] <- pars[[i]][pars[[i]] < ll] + 273.15
  }

  return(pars)
}


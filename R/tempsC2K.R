# Converts temperature parameters from C to K if values look like K (> cutoff, e.g., 270)

tempsC2K <- function(pars, cutoff) {

  for (i in grep('T_', names(pars))) {
    pars[[i]][pars[[i]] < cutoff] <- pars[[i]][pars[[i]] < cutoff] + 273.15
  }

  return(pars)
}


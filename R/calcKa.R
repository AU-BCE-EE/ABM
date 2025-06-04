# Calculate log Ka for speciation

calcKa <- function(pars) {

  # Caclulate numeic pKa from expressions and temp_K
  temp_K <- pars$temp_K
  pars$lkan <- as.numeric(lapply(pars$lka, function(x) eval(parse(text = x))))
  names(pars$lkan) <- names(pars$lka)

  return(pars)
}

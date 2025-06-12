# Returns COD per mol substrate

calcCOD <- function(form) {

  # If and only if first letter of form is lowercase, entire string is capitalized
  if(grepl('^[a-z]', form)) form <- toupper(form)
  # Read formula (function not vectorized)
  fc <- readFormula(form, elements = c('C', 'H', 'O', 'N'))
  # Calculate COD based on Rittmann and McCarty
  COD <- as.vector((2*fc['C'] + 0.5*fc['H'] - 1.5*fc['N'] - fc['O']) * molMass('O'))

  return(COD)
}

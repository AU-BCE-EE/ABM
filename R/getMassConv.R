# Get mass conversion factor to go from moles of component to mass COD, N, C, S, or total, in that order

getMassConv <- function(form) {

  # Remove p and m (+/-)
  form <- gsub('p$|m$', '', form)
  
  cod <- calcCOD(form)
  fn <- readFormula(form)
  
  if (cod > 0) {
    cf <- cod
  } else if ('N' %in% names(fn)) {
    cf <- molMass(form, elements = 'N')
  } else if ('C' %in% names(fn)) {
    cf <- molMass(form, elements = 'C')
  } else if ('S' %in% names(fn)) {
    cf <- molMass(form, elements = 'S')
  } else {
    cf <- molMass(form)
  }

  return(cf)

}

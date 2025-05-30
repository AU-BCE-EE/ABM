# Extract wash water volume

getWashWater <- function(pars) {
  
  if ('wash_water' %in% names(pars$var)) {
    wash_water <- pars$var[, 'wash_water']
  } else {
    wash_water <- 0 * pars$var[, 'slurry_mass']
  }

  return(wash_water)

}


getWashWater <- function(pars) {
  
  if ('wash_water' %in% names(pars$slurry_mass)) {
    wash_water <- pars$slurry_mass[, 'wash_water']
  } else {
    wash_water <- 0 * pars$slurry_mass[, 'slurry_mass']
  }

  return(wash_water)

}

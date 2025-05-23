
getRemVar <- function(pars, slurry_mass_approx) {
  # Determine slurry removal quantity in each time interval
  # Note final 0--alignment is a bit tricky
  if (slurry_mass_approx %in% c('late', 'mid')) {
    removals <- - c(0, 0, diff(pars$slurry_mass[-nrow(pars$slurry_mass), 'slurry_mass'])) > 0
  } else if (slurry_mass_approx == 'early') {
    removals <- - c(0, diff(pars$slurry_mass[, 'slurry_mass'])) > 0
  } 

  return(removals)
 
}

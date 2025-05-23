# Prepares slurry mass data for 

fixSlurryMass <- function(pars, days, slurry_mass_approx) {

  # Cannot have no slurry present because is used in all concentration calculations
  pars$slurry_mass[pars$slurry_mass[, 'slurry_mass'] == 0, 'slurry_mass'] <- 1E-10

  # Trim unused times
  pars$slurry_mass <- pars$slurry_mass[pars$slurry_mass$time <= days, ]

  # Check for sorted time
  if (is.unsorted(pars$slurry_mass$time)) {
    stop('Column `time` must be sorted when `slurry_mass` is time-dependent, but it is not: ',
         head(pars$slurry_mass$time))
  }
  
  # If simulation continues past slurry_mass data frame time, set following slurry production (addition) to zero
  # Extend last slurry mass all the way 
  if (pars$slurry_mass[nrow(pars$slurry_mass), 'time'] < days) {
    t_end <- days
    pars$slurry_mass <- rbind(pars$slurry_mass, pars$slurry_mass[nrow(pars$slurry_mass), ])
    pars$slurry_mass[nrow(pars$slurry_mass), 'time'] <- days
    # But make sure washing is not repeated!
    if (ncol(pars$slurry_mass) > 2) {
      pars$slurry_mass[nrow(pars$slurry_mass), 3:ncol(pars$slurry_mass)] <- 0
    }
  }

  if (slurry_mass_approx == 'mid') {
    # Get midpoint time
    ir <- which(- c(0, diff(pars$slurry_mass[, 'slurry_mass'])) > 0)
    tt <- (pars$slurry_mass[ir, 'time']  + pars$slurry_mass[ir - 1, 'time']) / 2
    # And copy slurry mass
    mm <- pars$slurry_mass[ir, 'slurry_mass']
    pars$slurry_mass <- rbind(pars$slurry_mass, data.frame(time = tt, slurry_mass = mm))
    pars$slurry_mass <- pars$slurry_mass[order(pars$slurry_mass$time), ]
  } 
  
  # Determine slurry removal quantity in each time interval
  # Note final 0--alignment is a bit tricky
  if (slurry_mass_approx %in% c('late', 'mid')) {
    removals <- - c(0, 0, diff(pars$slurry_mass[-nrow(pars$slurry_mass), 'slurry_mass'])) > 0
  } else if (slurry_mass_approx == 'early') {
    removals <- - c(0, diff(pars$slurry_mass[, 'slurry_mass'])) > 0
  } 

  pars$removals <- removals

  return(pars)

}

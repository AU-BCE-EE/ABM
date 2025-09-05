# Checks and prepares var par data, especially slurry_mass
# Applies approx_method to slurry_mass

fixVarDat <- function(pars, days) {

  if (is.null(pars$var)) {
    return(pars)
  }

  # Make sure at least slurry_mass is in pars$var data frame
  if (!'slurry_mass' %in% names(pars$var)) {
    stop('The pars var element is missing a slurry_mass column, which is required.')
  }

  # Cannot have no slurry present because is used in all concentration calculations
  pars$var[pars$var[, 'slurry_mass'] == 0, 'slurry_mass'] <- 1E-10

  # Trim unused times
  pars$var <- pars$var[pars$var$time <= days, ]

  # Check for sorted time
  if (is.unsorted(pars$var$time)) {
    stop('Column `time` must be sorted when time-variable parameters are used (pars), but it is not: ',
         head(pars$var$time))
  }

  # Check for duplicate times
  if (any(duplicated(pars$var$time))) {
    stop('Column `time` has duplicates, which is not allowed. See rows: ',
         paste0(which(duplicated(pars$var$time)), collape = ', '))
  }
  
  # If simulation continues past pars data frame time, extend last row all the way
  if (pars$var[nrow(pars$var), 'time'] < days) {
    t_end <- days
    pars$var <- rbind(pars$var, pars$var[nrow(pars$var), ])
    pars$var[nrow(pars$var), 'time'] <- days
    # But make sure washing is not repeated!
    if (ncol(pars$var) > 2) {
      pars$var[nrow(pars$var), 3:ncol(pars$var)] <- 0
    }
  }

  # For 'mid' option, other variables are copied from previous time
  if (pars$approx_method == 'mid') {
    # Get midpoint time
    ir <- which(- c(0, diff(pars$var[, 'slurry_mass'])) > 0)
    tt <- (pars$var[ir, 'time']  + pars$var[ir - 1, 'time']) / 2
    nr <- pars$var[ir, ]
    nr$time <- tt
    pars$var <- rbind(pars$var, nr)
    pars$var <- pars$var[order(pars$var$time), ]
  } 
  
  return(pars)

}

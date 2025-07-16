# Checks and prepares var par data, especially slurry_mass
# Applies approx_method to slurry_mass

fixVarDat <- function(var_pars, approx_method, days) {

  # Make sure at least slurry_mass is in var_pars$var data frame
  if (!'slurry_mass' %in% names(var_pars$var)) {
    stop('The var_pars var element is missing a slurry_mass column, which is required.')
  }

  # Cannot have no slurry present because is used in all concentration calculations
  var_pars$var[var_pars$var[, 'slurry_mass'] == 0, 'slurry_mass'] <- 1E-10

  # Trim unused times
  var_pars$var <- var_pars$var[var_pars$var$time <= days, ]

  # Check for sorted time
  if (is.unsorted(var_pars$var$time)) {
    stop('Column `time` must be sorted when time-variable parameters are used (var_pars), but it is not: ',
         head(var_pars$var$time))
  }
  
  # If simulation continues past var_pars data frame time, extend last row all the way
  if (var_pars$var[nrow(var_pars$var), 'time'] < days) {
    t_end <- days
    var_pars$var <- rbind(var_pars$var, var_pars$var[nrow(var_pars$var), ])
    var_pars$var[nrow(var_pars$var), 'time'] <- days
    # But make sure washing is not repeated!
    if (ncol(var_pars$var) > 2) {
      var_pars$var[nrow(var_pars$var), 3:ncol(var_pars$var)] <- 0
    }
  }

  # For 'mid' option, other variables are copied from previous time
  if (approx_method == 'mid') {
    # Get midpoint time
    ir <- which(- c(0, diff(var_pars$var[, 'slurry_mass'])) > 0)
    tt <- (var_pars$var[ir, 'time']  + var_pars$var[ir - 1, 'time']) / 2
    nr <- var_pars$var[ir, ]
    nr$time <- tt
    var_pars$var <- rbind(var_pars$var, nr)
    var_pars$var <- var_pars$var[order(var_pars$var$time), ]
  } 
  
  return(var_pars)

}

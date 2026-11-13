# Small internal helper functions and constants

# Order of conc_fresh elements as used in rates_cpp()
conc_fresh_names <- c('sulfide', 'urea', 'sulfate', 'TAN', 'starch', 'VFA', 'xa_dead', 'Cfat', 'CPs', 'CPf',
                      'RFd', 'iNDF', 'VSd', 'ash')

# Order of hydrolysis parameters (A and E) as used in rates_cpp()
hydrol_names <- c('xa_dead', 'starch', 'Cfat', 'CPs', 'CPf', 'RFd', 'VSd', 'urea')

# Order of COD_conv elements as used in rates_cpp()
COD_conv_names <- c('CH4', 'xa', 'RFd', 'iNDF', 'starch', 'Cfat', 'CP', 'VFA', 'S', 'VS', 'CO2_aer', 'CO2_sr',
                    'CO2_ureo', 'CP_N', 'C_xa', 'C_RFd', 'C_iNDF', 'C_starch', 'C_Cfat', 'C_CP', 'C_VFA',
                    'C_VSd', 'C_N_urea', 'C_xa_aer', 'C_xa_bac', 'C_xa_dead')

# Put named vector elements in the order expected by rates_cpp()
# Unnamed input is assumed to be in the correct order already
orderElements <- function(x, nms, what = '') {
  if (is.null(names(x))) {
    if (length(x) < length(nms)) stop('Parameter ', what, ' must have (named) elements ', paste(nms, collapse = ', '))
    return(x)
  }
  if (any(miss <- !nms %in% names(x))) {
    stop('Parameter ', what, ' is missing element(s): ', paste(nms[miss], collapse = ', '))
  }
  # Keep any extra elements at the end
  x <- unlist(c(x[nms], x[!names(x) %in% nms]))
  return(x)
}

# Convert time-variable input (numeric constant or data frame with time in column 1 and value in column 2)
# to standard data frame with columns time and y (or return constant unchanged)
timeDF <- function(dat) {
  if (is.data.frame(dat)) {
    return(data.frame(time = dat[[1]], y = dat[[2]]))
  }
  if (!is.numeric(dat)) stop('Time-variable input must be numeric or a data frame.')
  return(dat)
}

# Keep only some rows (output times) of simulation output
# Effluent (columns ending in _eff) in dropped rows is added to the next kept row, so cumulative effluent is correct
filterTimes <- function(dat, keep) {
  if (all(keep)) return(dat)
  eff_cols <- grep('_eff$', names(dat), value = TRUE)
  idx_keep <- which(keep)
  if (length(eff_cols) > 0 && length(idx_keep) > 0) {
    # Row index of first kept row at or after each row
    target <- idx_keep[findInterval(seq_along(keep), idx_keep, left.open = TRUE) + 1]
    ok <- !is.na(target)
    eff <- as.matrix(dat[ok, eff_cols, drop = FALSE])
    eff[is.na(eff)] <- 0
    eff_sum <- rowsum(eff, target[ok])
    dat[as.integer(rownames(eff_sum)), eff_cols] <- eff_sum
  }
  return(dat[keep, , drop = FALSE])
}

# Average rate from cumulative values, e.g., emission rate from cumulative emission
# For duplicated times (e.g., at emptying) the previous rate is used
avgRate <- function(cum, time) {
  dt <- c(1, diff(time))
  r <- c(0, diff(cum)) / dt
  for (i in which(dt == 0)) r[i] <- r[i - 1]
  return(r)
}

# Run deSolve::lsoda() for one interval
# compiled = TRUE: derivatives from the compiled function rates_derivs() (see rates_cpp.cpp), which deSolve calls
# directly, so no R code runs during integration. Parameters are passed to C++ once with rates_set_pars().
# compiled = FALSE: rates_cpp() is called as an R function (slower, same results, for testing and debugging)
# Set options(ABM.compiled = FALSE) to use the R function in abm()
runLsoda <- function(y, times, pars, compiled = getOption('ABM.compiled', TRUE)) {
  if (compiled) {
    dims <- rates_set_pars(pars)
    if (length(y) != dims[1]) stop('Length of y (', length(y), ') does not match number of state variables (', dims[1], ')')
    out <- deSolve::lsoda(y = y, times = times, func = 'rates_derivs', dllname = 'ABM', initfunc = NULL,
                          parms = NULL, nout = dims[2], outnames = pars$out_names, rtol = 1E-5, atol = 1E-5)
  } else {
    out <- deSolve::lsoda(y = y, times = times, func = rates_cpp, parms = pars, rtol = 1E-5, atol = 1E-5)
  }
  return(out)
}

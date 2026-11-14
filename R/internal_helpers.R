# Small internal helper functions and constants

# Order of conc_fresh elements as used in rates_cpp()
# VFA_C (g C/kg) is optional in conc_fresh, default is VFA as acetate (set in abm())
conc_fresh_names <- c('sulfide', 'urea', 'sulfate', 'TAN', 'starch', 'VFA', 'xa_dead', 'Cfat', 'CPs', 'CPf',
                      'RFd', 'iNDF', 'VSd', 'ash', 'VFA_C')

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

# Put a fixed (not data frame) value from add_pars into a time-variable parameter given as a data frame
# (e.g., conc_fresh, xa_fresh, temp_C). The data frame structure is kept: the given columns become constant
# over time, the time column and all other columns are unchanged.
# val: named values set those columns; an unnamed value sets the column named like the parameter (nm) if there
#   is one (e.g., temp_C, pH, slurry_mass), otherwise all value columns (several unnamed values: by position)
# grp: microbial group parameter, where names 'default' (groups without own value) and 'all' (every group) are keywords
fixedIntoDF <- function(df, val, nm, grp = FALSE, warn = TRUE) {
  vcols <- setdiff(names(df), 'time')
  v <- unlist(val)
  v_nms <- names(v)
  if (is.null(v_nms) || all(v_nms == '')) {
    target <- if (nm %in% vcols) nm else vcols
    if (length(v) == 1) {
      df[target] <- v[[1]]
    } else if (length(v) == length(target)) {
      df[target] <- as.list(unname(v))
    } else {
      stop('`add_pars` element ', nm, ' is unnamed and has length ', length(v), ', but the time-variable parameter has ',
           length(target), ' value column(s) (', paste(target, collapse = ', '), '). Give a single value, one value per column, or named values.')
    }
    return(df)
  }
  if (any(v_nms == '')) stop('`add_pars` element ', nm, ' has some unnamed values. Name all values or none.')
  if ('time' %in% v_nms) stop('`add_pars` element ', nm, ': the time column of a time-variable parameter cannot be set.')
  kw <- if (grp) c('all', 'default') else character()
  named <- setdiff(v_nms, kw)
  if ('default' %in% v_nms && grp) df[vcols] <- v[['default']]
  if (length(named) > 0) df[named] <- as.list(v[named])
  if ('all' %in% v_nms && grp) {
    if (warn && length(named) > 0) {
      warning('add_pars element ', nm, ' has both an "all" value and group-specific values (', paste(named, collapse = ', '),
              '). "all" sets every group, so the group-specific values are ignored.\n  Use "default" instead of "all" to give a value for the groups without their own value.')
    }
    df[setdiff(names(df), 'time')] <- v[['all']]
  }
  df
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

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

# Creates an interpolation table for time-variable inputs that is used directly in rates_cpp()
# (replaces R interpolation functions, which are slow to call from C++)
# Interpolation in rates_cpp() is identical to approxfun(..., rule = 2, f = 0)

# dat = numeric constant (vector or list, possibly named) or data frame with column `time`
# cols = names of columns (or elements for constant input) to return, in this order
# approx_method = 'linear' or 'constant'

makeInterpTab <- function(dat, cols = NULL, approx_method = 'linear') {

  method <- match(approx_method, c('linear', 'constant')) - 1L
  if (is.na(method)) {
    stop('approx_method must be "linear" or "constant", but is "', approx_method, '".')
  }

  if (is.data.frame(dat)) {
    if (!'time' %in% names(dat)) stop('Time-variable input must have a column named `time`.')
    if (is.null(cols)) cols <- names(dat)[names(dat) != 'time']
    if (any(miss <- !cols %in% names(dat))) {
      stop('Missing columns in time-variable input: ', paste(cols[miss], collapse = ', '))
    }
    dat <- dat[!is.na(dat$time), , drop = FALSE]
    dat <- dat[order(dat$time), , drop = FALSE]
    y <- as.matrix(dat[, cols, drop = FALSE])
    storage.mode(y) <- 'double'
    time <- as.numeric(dat$time)
    # Duplicate times are collapsed to mean values (as in approxfun())
    if (anyDuplicated(time)) {
      warning('Duplicated times in time-variable input; collapsing to unique times (mean values).')
      y <- apply(y, 2, function(x) tapply(x, time, mean))
      y <- matrix(y, ncol = length(cols))
      time <- sort(unique(time))
    }
    if (anyNA(y)) stop('Missing values (NA) in time-variable input: ', paste(cols, collapse = ', '))
  } else {
    v <- unlist(dat)
    if (!is.null(cols)) {
      if (any(miss <- !cols %in% names(v))) stop('Missing elements: ', paste(cols[miss], collapse = ', '))
      v <- v[cols]
    }
    y <- matrix(as.numeric(v), nrow = 1)
    time <- 0
  }

  colnames(y) <- cols
  return(list(time = time, y = y, method = method))

}

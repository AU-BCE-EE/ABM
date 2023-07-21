# Creates a function for time-variable inputs

makeTimeFunc <- function(dat, x = 1, y = 2, approx_method = 'linear') {

  if (is.data.frame(dat)) {
    x <- dat[, x]
    y <- dat[, y]
    func <- approxfun(x, y, method = approx_method, rule = 2)
  } else if (is.numeric(dat)) {
    func <- function(x) return(dat)
  } else {
    stop('Input to makeTimeFunc must be numeric or data frame.')
  }

  return(func)

}


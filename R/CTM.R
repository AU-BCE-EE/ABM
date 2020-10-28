# Temperature model

CTM <- function(tt, topt, tmin, tmax, yopt = 1) {

  if (!all.equal(length(topt), length(tmin), length(tmax), length(yopt))) {
    stop('Length of topt and following arguments be be identical.')
  }

  if (length(tt) > 1) {
    tt <- tt[1]
    warning('Multiple temperature (tt) values given but only first will be used')
  }

  if (length(topt) > 1) {
    y <- NULL
    for (i in 1:length(topt)) {
      y <- c(y, CTM(tt, topt[i], tmin[i], tmax[i], yopt[i]))
    }

    names(y) <- names(topt)
    return(y)

  }

  # Scale (really only needed for flip)
  tt <- tt - topt
  tmax <- tmax - topt
  tmin <- tmin - topt
  topt <- topt - topt
  
  flip <- FALSE

  if (topt - tmin < (tmax - tmin)/2) {
    flip <- TRUE
    tt <- - tt
    tmaxo <- tmax
    tmino <- tmin
    tmax <- - tmino
    tmin <- - tmaxo
  }
  
  y <- yopt * ((tt - tmax) * (tt - tmin)^2) / 
          ((topt - tmin) * ((topt - tmin) * (tt - topt) - 
                              (topt - tmax) * (topt + tmin - 2*tt))
           )
  
  if (flip) {
    tt <- - tt
    tmin <- tmino
    tmax <- tmaxo
  }
  
  # Fix values outside (or close to) limits
  y[tt <= tmin | tt >= tmax] <- 0
  y[y < 0] <- 0

  return(y)
}



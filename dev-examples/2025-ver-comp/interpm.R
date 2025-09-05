# Replace missing values in multiple data frame columns with interpolated values

interpm <- function(dat, x, ys, by = NA, ...) {

  if (is.na(by)) {
    for (i in ys) {
      rout <- which(is.na(dat[[i]])) 
      if (length(rout) > 0) {
        dat[[i]][rout] <- approx(dat[[x]][-rout], dat[[i]][-rout], xout = dat[[x]][rout], ...)$y 
      } 
    }
  } else {
    for (i in ys) {
      for (j in unique(dat[[by]])) {
        gr <- dat[[by]] == j
        rout <- is.na(dat[[i]])
        if (length(rout) > 0) {
          dat[[i]][gr & rout] <- approx(dat[[x]][gr & !rout], dat[[i]][gr & !rout], xout = dat[[x]][gr & rout], ...)$y 
        } 
      }
    }

  }

  return(dat)

}

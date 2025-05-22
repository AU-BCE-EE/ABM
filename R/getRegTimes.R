

getRegTimes <- function(pars, days) {
  
  # If empty interval is set to 0 or NA the storage is never emptied. 
  empty_int <- pars$empty_int
  if(empty_int == 0 || is.na(empty_int)) {
    empty_int <- days + 1
  }
  
  # Figure out time intervals for loop
  if (!is.na(pars$wash_int) && pars$wash_water > 0) {  
    wash_int <- pars$wash_int
    rest_d <- pars$rest_d
  } else {
    wash_int <- Inf
    rest_d <- 0
  }
  wash_rest_int <- wash_int + rest_d

  # Continue sorting out intervals
  i <- 0
  t_int <- 0
  t_nowash <- 0
  wash <- FALSE

  # Continute . . .
  # Each interval is either 1) the fixed empty_int or if time between washings would be exceeded, 2) time to get to a washing event
  while (sum(t_int, wash * rest_d) < days) {
    i <- i + 1
    t_int[i] <- min(wash_int - t_nowash, empty_int, days - sum(t_int, wash * rest_d))
    if (t_int[i] == wash_int - t_nowash) {
      wash[i] <- TRUE
      t_nowash <- 0
    } else {
      wash[i] <- FALSE
      t_nowash <- t_nowash + t_int[i]
    }
  }

  # Number of empty or wash intervals
  n_int <- length(t_int)

  return(list(n_int = n_int,
              empty_int = empty_int,
              wash_int = wash_int,
              rest_d = rest_d,
              t_int = t_int,
              wash = wash))

}

# Internal function for repeated abm() calls
# Argument list should exactly match abm()

abmStartup <- function(days,
                      delta_t,
                      times,
                      pars,
                      startup,
                      starting,
                      value,
                      warn) {

  value.orig <- value
  value <- 'ts'

  for (i in 1:(startup + 1)) {
    if (i > startup) {
      cat('and final run')
      cat('\n')
    } else {
      cat(paste0(i, 'x -> '))
    }

    if (i > startup) {
      value <- value.orig
    }

    print(i)
    if (i > 1) {
      starting <- out
    browser()
    } 

    
    # Call abm() with arguments given in outside call except for startup and value
    out <- abm(days = days,
               delta_t = delta_t,
               times = times,
               pars = pars,
               startup = 0,
               starting = starting,
               value = value,
               warn = warn)

  }

  return(out)

}

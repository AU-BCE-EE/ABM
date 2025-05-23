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

    if (i > 1) {
      starting <- out
    } 
    
    # Call abm() with arguments given in outside call except for startup and value
    out <- abm(days = days,
               delta_t = delta_t,
               times = times,
               pars = pars,
               startup = 0,
               starting = starting,
               value = 'ts',
               warn = warn)
   
    if (i <= startup) {
      # Pull starting *concentrations* (inlcuding xa) from previous sim
      tso <- out

      # Names need to deal with possible data frame for conc_fresh
      cf_names <- names(man_pars$conc_fresh)
      cf_names <- cf_names[!grepl('^time', paste0(cf_names, '_conc'))]

      init_pars$conc_init <- unlist(tso[nrow(tso), paste0(cf_names, '_conc')])
      names(init_pars$conc_init) <- cf_names

      grp_pars$xa_init <- unlist(tso[nrow(tso), paste0(grp_pars$grps, '_conc')])
      names(grp_pars$xa_init) <- grp_pars$grps
    }
 
  }

  return(out)

}

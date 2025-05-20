# Internal function for repeated abm() calls
# Argument list should exactly match abm()

abmRepeat <- function(days,
                      delta_t,
                      times,
                      wthr_pars,
                      evap_pars,
                      mng_pars,
                      man_pars,
                      init_pars,
                      grp_pars,
                      mic_pars,
                      chem_pars,
                      arrh_pars,
                      anim_pars,
                      resp,
                      pH_inhib_overrule,
                      add_pars,
                      pars,
                      startup,
                      starting,
                      approx_method,
                      par_key,
                      value,
                      rates_calc,
                      warn)

  # Check for conc_init pars in add_pars--this is not compatible with startup
  if (any(grepl('conc_init', names(add_pars)))) {
    stop('Simulation has a startup period (startup > 0) and initial concentrations in add_pars.\n  These two options do not work together--see issue #57')
  }

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

    # Call abm() with arguments given in outside call except for startup and value
    out <- abm(days = days,                     days,              
               delta_t = delta_t,               delta_t,                          
               times = times,                   times,
               wthr_pars = wthr_pars,           wthr_pars,
               evap_pars = evap_pars,           evap_pars,
               mng_pars = mng_pars,             mng_pars,
               man_pars = man_pars,             man_pars,
               init_pars = init_pars,           init_pars,
               grp_pars = grp_pars,             grp_pars,
               mic_pars = mic_pars,             mic_pars,
               chem_pars = chem_pars,           chem_pars,
               arrh_pars = arrh_pars,           arrh_pars,
               anim_pars,
               add_pars = add_pars,             resp,
               pars = pars,                     pH_inhib_overrule,
                                                add_pars,
               startup = 0,                     pars,
               starting = starting,             startup,
               approx_method = approx_method,   starting,
               par_key = par_key,               approx_method,
               value = value,                   par_key,
               warn = warn)                     value,
                                                rates_calc,
                                                warn)
                                                
    
 
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

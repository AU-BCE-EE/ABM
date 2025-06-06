makeInitState <- function(pars) {
 
  if(is.data.frame(pars$slurry_mass)){
    # If missing slurry mass at time 0, set to earliest slurry mass value
    if (pars$slurry_mass[1, 'time'] > 0) {
      pars$slurry_mass <- rbind(c(0, pars$slurry_mass$slurry_mass[1]), pars$slurry_mass)
    }
    slurry_mass_init <- pars$slurry_mass[1, 'slurry_mass']
  } else {
    slurry_mass_init <- pars$slurry_mass
  }
  
  if (slurry_mass_init == 0) {
    slurry_mass_init <- 1E-10
  }

  if (!is.null(pars$kla)) {
    emis <- pars$kla * 0
    names(emis) <- paste0(names(emis), '_emis_cum')
  } else {
    emis <- NULL
  }
  
  y <- c(pars$xa_init * slurry_mass_init,                       # Multiple microbial groups
         pars$sub_init[pars$subs] * slurry_mass_init,           # Multiple particulate substrates 
         pars$conc_init * slurry_mass_init,                     # VFA and conservative solutes
         slurry_mass = slurry_mass_init, 
         CH4_emis_cum = 0, 
         CO2_emis_cum = 0, 
         emis,
         slurry_load = 0,
         COD_load = 0)
 
  return(y)

}

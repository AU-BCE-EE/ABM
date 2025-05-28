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
  
  y <- c(pars$xa_init * slurry_mass_init,                        # Multiple microbial groups
         pars$conc_init[pars$subs] * slurry_mass_init,           # Multiple particulate substrates 
         VFA = pars$conc_init[['VFA']] * slurry_mass_init,       # VFA
         slurry_mass = slurry_mass_init, 
         CH4_emis_cum = 0, 
         CO2_emis_cum = 0, 
         slurry_load = 0,
         COD_load = 0)
 
  return(y)

}

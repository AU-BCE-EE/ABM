makeInitState <- function(pars, 
                          starting) {
 
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
  
  y <- c(xa = pars$xa_init * slurry_mass_init,
         slurry_mass = slurry_mass_init, 
         xa_aer = pars$conc_init[['xa_aer']] * slurry_mass_init,
         xa_bac = pars$conc_init[['xa_bac']] * slurry_mass_init,
         xa_dead = pars$conc_init[['xa_dead']] * slurry_mass_init, 
         RFd = pars$conc_init[['RFd']] * slurry_mass_init,
         iNDF = pars$conc_init[['iNDF']] * slurry_mass_init,
         ash = pars$conc_init[['ash']] * slurry_mass_init,
         VSd = pars$conc_init[['VSd']] * slurry_mass_init,
         starch = pars$conc_init[['starch']] * slurry_mass_init,
         CPs = pars$conc_init[['CPs']] * slurry_mass_init,
         CPf = pars$conc_init[['CPf']] * slurry_mass_init,
         Cfat = pars$conc_init[['Cfat']] * slurry_mass_init,
         VFA = pars$conc_init[['VFA']] * slurry_mass_init, 
         urea = pars$conc_init[['urea']] * slurry_mass_init, 
         TAN = pars$conc_init[['TAN']] * slurry_mass_init, 
         sulfate = pars$conc_init[['sulfate']] * slurry_mass_init, 
         sulfide = pars$conc_init[['sulfide']] * slurry_mass_init, 
         NH3_emis_cum = 0, 
         N2O_emis_cum = 0, 
         CH4_emis_cum = 0, 
         CO2_emis_cum = 0, 
         COD_conv_cum = 0, 
         COD_conv_cum_meth = 0, 
         COD_conv_cum_respir = 0, 
         COD_conv_cum_sr = 0,
         COD_load_cum = 0,
         C_load_cum = 0,
         N_load_cum = 0,
         slurry_load_cum = 0)
  
  if (!is.null(starting) & is.data.frame(starting)) {
    start.vars <- c('slurry_mass', 'xa_aer', 'xa_bac', 'xa_dead', 'iNDF', 'ash', 'RFd', 'VSd', 'starch', 'CPs', 'CPf', 'Cfat', 'VFA', 'urea', 'TAN', 'sulfate', 'sulfide')
    y[start.vars]  <- starting[nrow(starting), start.vars]
  }  

  return(y)

}

rates <- function(t, 
		              y, 
		              parms, 
		              temp_C_fun = temp_C_fun, 
		              pH_fun = pH_fun) {
    
  # Short name for parms to make indexing in code below simpler 
  p <- parms

  # Time-dependent values
  pH <- pH_fun(t + p$t_run)
  temp_C <- temp_C_fun(t + p$t_run)
  temp_K <- temp_C + 273.15

  # Hydrolysis rate for now
  alpha <- CTM_cpp(temp_K, p$T_opt_hyd, p$T_min_hyd, p$T_max_hyd, p$hydrol_opt)

  # Create vectors with derivative components, all with order of y elements
  qhat <- rut <- consump <- growth <- inflow <- death <- hydrol <- emis <- 0 * y

  # Microbial substrate utilization rate (vectorized calculation)
  qhat[p$i_mic] <- CTM_cpp(temp_K, p$T_opt, p$T_min, p$T_max, p$qhat_opt)

  # VFA consumption rates (g/d) and growth
  rut[p$i_meth] <- (qhat[p$i_meth] * y['VFA'] * y[p$i_meth]) / (p$ks[p$i_meth] * y['slurry_mass'] + y['VFA'])
  rut[p$i_sr] <- qhat[p$i_sr] * 0
  growth[p$i_mic] <- p$yield * rut[p$i_mic]
  consump['VFA'] <- - sum(rut[p$i_meth])

  # Inflow from slurry addition
  inflow[p$i_mic] <- p$xa_fresh * p$slurry_prod_rate
  inflow[c('VFA', 'VSd')] <- p$conc_fresh[c('VFA', 'VSd')] * p$slurry_prod_rate
  inflow[c('slurry_mass', 'slurry_load')] <- p$slurry_prod_rate
  inflow['COD_load'] <- sum(p$conc_fresh[c('VSd', 'VFA')], p$xa_fresh) * p$slurry_prod_rate

  # Death of microbes
  death[p$i_mic] <- - p$dd_rate * y[p$i_mic]
  death['VSd'] <- - sum(death[p$i_mic])

  # Hydrolysis
  hydrol[c('VSd', 'VFA')] <- c(-1, 1) * alpha * y['VSd']

  # Emission
  p$COD_conv['CO2_meth'] <- 5
  emis[c('CH4_emis_cum', 'CO2_emis_cum')] <- sum(rut[p$i_meth]) / p$COD_conv[c('CH4', 'CO2_meth')]

  # Add vectors to get derivatives
  # All elements in g/d as COD except slurry_mass (kg/d as fresh slurry mass)
  ders <- inflow + growth + consump + death + hydrol + emis
 
  return(list(ders, c(CH4_emis_rate = emis[['CH4_emis_cum']], temp_C = temp_C, pH = pH)))

}

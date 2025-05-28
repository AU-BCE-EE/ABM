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
  names(alpha) <- p$subs

  # Create vectors with derivative components, all with order of y elements
  qhat <- rut <- consump <- growth <- inflow <- death <- hydrol <- emis <- 0 * y

  # Microbial substrate utilization rate (vectorized calculation)
  qhat[p$grps] <- CTM_cpp(temp_K, p$T_opt, p$T_min, p$T_max, p$qhat_opt)

  # VFA consumption rates (g/d) and growth
  rut[p$meths] <- (qhat[p$meths] * y['VFA'] * y[p$meths]) / (p$ks[p$meths] * y['slurry_mass'] + y['VFA'])
  rut[p$srs] <- qhat[p$srs] * 0
  growth[p$grps] <- p$yield * rut[p$grps]
  consump['VFA'] <- - sum(rut[p$meths])

  # Inflow from slurry addition
  inflow[p$grps] <- p$xa_fresh * p$slurry_prod_rate
  inflow[p$subs] <- p$conc_fresh[p$subs] * p$slurry_prod_rate
  inflow['VFA'] <- p$conc_fresh['VFA'] * p$slurry_prod_rate
  inflow[c('slurry_mass', 'slurry_load')] <- p$slurry_prod_rate
  inflow['COD_load'] <- sum(p$conc_fresh[c('VSd', 'VFA')], p$xa_fresh) * p$slurry_prod_rate

  # Death of microbes
  death[p$grps] <- - p$dd_rate * y[p$grps]
  death['VSd'] <- sum(death[p$grps])

  # Hydrolysis
  hydrol[p$subs] <- -alpha * y[p$subs]
  hydrol['VFA'] <- - sum(hydrol[p$subs])

  # Emission
  p$COD_conv['CO2_meth'] <- 5
  emis[c('CH4_emis_cum', 'CO2_emis_cum')] <- sum(rut[p$meths]) / p$COD_conv[c('CH4', 'CO2_meth')]

  # Add vectors to get derivatives
  # All elements in g/d as COD except slurry_mass (kg/d as fresh slurry mass)
  ders <- inflow + growth + consump + death + hydrol + emis
 
  return(list(ders, c(CH4_emis_rate = emis[['CH4_emis_cum']], temp_C = temp_C, pH = pH)))

}

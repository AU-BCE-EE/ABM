rates <- function(t, 
		              y, 
		              parms) {
    
  # Short name for parms to make indexing in code below simpler 
  p <- parms

  # Calculation speciation, used for emission and inhibition
  p <- calcSpec(p, y)

  # Determine inhibition reductions
  p <- calcInhib(p, y)

  # Initialize vectors with derivative components, all with same order of y elements
  alpha <- qhat <- rut <- consump <- growth <- inflow <- death <- hydrol <- volat <- emis <- 0 * y

  # Hydrolysis rate
  alpha[p$subs] <- CTM_cpp(p$temp_K, p$T_opt_hyd, p$T_min_hyd, p$T_max_hyd, p$hydrol_opt)

  # Microbial substrate utilization rate (vectorized calculation)
  qhat[p$grps] <- CTM_cpp(p$temp_K, p$T_opt, p$T_min, p$T_max, p$qhat_opt)

  # VFA consumption rates (g/d) and growth
  # Methanogens
  rut[p$meths] <- p$ired[p$meths] * (qhat[p$meths] * y['VFA'] * y[p$meths]) / (p$ks[p$meths] * y['slurry_mass'] + y['VFA'])
  # Sulfate reducers
  rut[p$srs] <- p$ired[p$srs] * qhat[p$srs] * 0
  
  # VFA consumption is sum of all rut terms
  consump['VFA'] <- - sum(rut)
  
  # Growth rate of all groups
  growth[p$grps] <- p$yield[p$grps] * rut[p$grps]

  # Inflow from slurry addition
  # First only concentrations are set, and multiplied by inflow in last line
  inflow[p$grps] <- p$xa_fresh
  inflow[p$subs] <- p$sub_fresh[p$subs]
  inflow[p$sols] <- p$conc_fresh[p$sols]
  inflow[c('slurry_mass', 'slurry_load')] <- 1
  inflow['COD_load'] <- sum(inflow[p$grps], inflow[p$subs] * p$stoich['VFA', p$subs], inflow['VFA'])
  inflow <- inflow * p$slurry_prod_rate

  # Death of microbes
  death[p$grps] <- - p$dd_rate * y[p$grps]
  # NTS: we need an input parameter setting the sink for dd to a certain substrate
  # NTS: could be one just for the purpose, like xa_dead
  death[p$subs[1]] <- - sum(death[p$grps])

  # Hydrolysis
  hydrol[p$subs] <- - alpha[p$subs] * y[p$subs]
  # Production of arbitrary components based on specified stoichiometry (can omit components)
  hydrol[rownames(p$stoich)] <- - p$stoich %*% hydrol[colnames(p$stoich)]
  

  # Volatilization
  volat <- calcVolat(p, volat)
  
  # Emission of CH4 and CO2
  p$COD_conv['CO2_meth'] <- 5
  emis[c('CH4_emis_cum', 'CO2_emis_cum')] <- (sum(rut[p$meths]) - sum(growth[p$meths])) / 
                                              p$COD_conv[c('CH4', 'CO2_meth')]
  # Add vectors to get derivatives
  # All elements in g/d as COD except 
  #   * slurry_mass (kg/d as fresh slurry mass)
  #   * CH4 (g/d as CH4 or C?)
  #   * solutes other than VFA (g/d)
  ders <- inflow + growth + consump + death + hydrol + volat + emis

  return(list(ders, c(CH4_emis_rate = emis[['CH4_emis_cum']], temp_C = p$temp_C, pH = p$pH)))

}

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
  rut <- respir <- consump <- growth <- inflow <- death <- hydrol <- volat <- meth <- 0 * y

  # Other (temperature-dependent) derivative vectors, brought in with pars
  alpha <- p$alpha
  qhat <- p$qhat

  # VFA consumption rates (g/d) and growth
  # Rate of substrate utilization
  rut[p$meths] <- p$ired[p$meths] * qhat[p$meths] * y[p$meths] * y['CH3COOH'] / 
                  (p$ksv[p$meths] * y['slurry_mass'] + y['CH3COOH'])
               
  # Sulfate and H2S
  # NTS: missing S:COD conversion factor!
  if (!p$sromit) {
    rut[p$srs] <- p$ired[p$srs] * qhat[p$srs] * y[p$srs] * 
                  y['CH3COOH'] / (p$ksv[p$srs] * y['slurry_mass'] + y['CH3COOH']) *
                  y['SO4m2'] / (p$kss[p$srs] * y['slurry_mass'] + y['SO4m2'])
    consump['SO4m2'] <- - sum(rut[p$srs] * (1 - p$yield[p$srs]))
    if ('H2S' %in% names(y)) {
      consump['H2S'] <- - consump['SO4m2']
    }
  }
  
  # VFA consumption is sum of all rut terms
  consump['CH3COOH'] <- - sum(rut)
 
  # Growth rate of all groups
  growth[p$grps] <- p$yield[p$grps] * rut[p$grps]

  # Respiration only consumes VFA
  if (!is.null(p$O2kl)) {
    respir['CH3COOH'] <- - 0.01 * p$O2kl * p$area
    if ('CO2' %in% names(y)) {
      respir['CO2'] <- - 12.01 / 32. * respir['CH3COOH']
    }
  }

  # Inflow from slurry addition
  # First only concentrations are set, and multiplied by inflow in last line
  inflow[p$grps] <- p$xa_fresh
  inflow[p$subs] <- p$sub_fresh[p$subs]
  inflow[p$sols] <- p$conc_fresh[p$sols]
  inflow[c('slurry_mass', 'slurry_load')] <- 1
  inflow['COD_load'] <- sum(inflow[p$grps], inflow[p$subs] * p$stoich['CH3COOH', p$subs], inflow['CH3COOH'])
  inflow <- inflow * p$slurry_prod_rate

  # Death of microbes
  death[p$grps] <- - p$dd_rate[p$grps] * y[p$grps]
  # NTS: we need an input parameter setting the sink for dd to a certain substrate
  # NTS: could be one just for the purpose, like xa_dead
  death[p$subs[1]] <- - sum(death[p$grps])

  # Hydrolysis
  hydrol[p$subs] <- - alpha[p$subs] * y[p$subs]
  # Production of arbitrary components based on specified stoichiometry (can omit components)
  hydrol[rownames(p$stoich)] <- - p$stoich %*% hydrol[colnames(p$stoich)]
  
  # Volatilization (can include CO2)
  volat <- calcVolat(p, volat)
  
  # Methanogenesis
  # All CH4 emitted
  meth['CH4_emis_cum'] <- (sum(rut[p$meths]) - sum(growth[p$meths])) / p$COD_conv['CH4']
  # CO2 goes into dissolved pool
  # Note fixed coefficient 1:1 CO2:CH4 on a C basis, because all methane is always from acetic acid
  if ('CO2' %in% names(y)) {
    meth['CO2'] <- meth['CO2'] + meth['CH4_emis_cum'] * 1 
  }
  # Add vectors to get derivatives
  # All elements in g/d as COD except 
  #   * slurry_mass (kg/d as fresh slurry mass)
  #   * CH4 (g/d as CH4 or C?)
  #   * solutes other than VFA (...)
  ders <- inflow + growth + respir + consump + death + hydrol + volat + meth

  return(list(ders, c(CH4_emis_rate = meth[['CH4_emis_cum']], temp_C = p$temp_C, pH = p$pH)))

}

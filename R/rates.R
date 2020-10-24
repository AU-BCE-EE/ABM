rates <-
function(t, y, parms, temp_C_fun = temp_C_fun, pH_fun = pH_fun, SO4_fun = SO4_fun) {

  y[y < 0] <- 0

  ks_SO4 <- parms$ks_SO4
  ki_H2S_meth <- parms$ki_H2S_meth
  ki_H2S_sr <- parms$ki_H2S_sr
  alpha_opt <- parms$alpha_opt
  alpha_T_opt <- parms$alpha_T_opt
  alpha_T_min <- parms$alpha_T_min
  alpha_T_max <- parms$alpha_T_max
  
  slurry_prod_rate <- parms$slurry_prod_rate               
  max_slurry_mass <- parms$ max_slurry_mass
  resid_frac <- parms$resid_frac
  area <- parms$area
  temp_C <- parms$temp_C
                  
  conc_fresh <- parms$conc_fresh
  
  yield <- parms$yield
  xa_fresh <- parms$xa_fresh
  xa_init <- parms$xa_init
  decay_rate <- parms$decay_rate 
  ks_coefficient <- parms$ks_coefficient
  resid_enrich <- parms$resid_enrich
  qhat_opt <- parms$qhat_opt  
  T_opt <- parms$T_opt
  T_min <- parms$T_min
  T_max <- parms$T_max
  ki_NH3_min <- parms$ki_NH3_min
  ki_NH3_max <- parms$ki_NH3_max
  ki_NH4_min <- parms$ki_NH4_min
  ki_NH4_max <- parms$ki_NH4_max
  pH_upr <- parms$pH_upr
  pH_lwr <- parms$pH_lwr
  
  COD_conv <- parms$COD_conv
  kl<- parms$kl  
  
  t_run <- parms$t_run
  
  # Hard-wired settings
  temp_standard <- 298
  temp_zero <- 273

  # Temp in K
  temp_C <- temp_C_fun(t + t_run)
  temp_K <- temp_C + 273.15

  # Find methanogens
  i_meth <- which(grepl('^[mp]', names(qhat_opt)))
  i_sr <- which(grepl('^sr', names(qhat_opt)))
  n_mic <- length(qhat_opt)

  # Extract state variable values from y argument
  # Slurry mass (kg)
  slurry_mass <- y[['slurry_mass']]
  # Particulate substrate mass Sp
  Sp <- y[['Sp']]
  # VFA mass
  VFA <- y[['VFA']]
  # Xa mass (g)
  xa <- y[1:n_mic]
  # Sulfate mass (g)
  sulfate <-y[['sulfate']]
  # Sulfide
  sulfide <- y[['sulfide']]
  # Cumulative methane production
  CH4_emis_cum <- y[['CH4_emis_cum']]

  # pH and SO4-2
  conc_fresh_SO4 <- SO4_fun(t + t_run)
  conc_SO4 <- sulfate/slurry_mass
  if (is.numeric(parms$pH) | is.data.frame(parms$pH)) {
    pH <- pH_fun(t + t_run)
  } else if (parms$pH == 'calc') {
    pH <- H2SO4_titrat(conc_SO4)
  } else {
    stop('pH problem (xi342)')
  }

  # Hard-wired equilibrium constants
  log_ka <- c(NH3 = - 0.0905 - 2729.31/temp_K, 
            H2S = -7.051 + exp(0.029 * (temp_K - temp_standard)))

  # Hard-wired Henry's law constants
  kH_oxygen <- 0.0013*exp(1700*((1/temp_K)-(1/temp_standard)))*32*1000 # g/(m3*abm)
  
  # Derived parameters
  # Hydrolysis rate
  alpha <- alpha_opt * ((temp_K - alpha_T_max) * (temp_K - alpha_T_min)^2) / 
             ((alpha_T_opt - alpha_T_min) * ((alpha_T_opt - alpha_T_min) * (temp_K - alpha_T_opt) - 
                                            (alpha_T_opt - alpha_T_max) * 
                                            (alpha_T_opt + alpha_T_min - 2*temp_K))
             )

  if (alpha < 0) alpha <- 0
  
  # Microbial substrate utilization rate (vectorized calculation)
  qhat <- qhat_opt * ((temp_K - T_max) * (temp_K - T_min)^2) / 
            ((T_opt - T_min) * ((T_opt - T_min) * (temp_K - T_opt) - 
                                (T_opt - T_max) * (T_opt + T_min - 2*temp_K))
            )
  # Set to zero outside bounds
  qhat[temp_K < T_min | temp_K > T_max] <- 0
  
  # Ks temperature dependence
  ks <- ks_coefficient*(0.8157*exp(-0.063*temp_C))#(71.3*exp(-0.175*temp_C))

  # Chemical speciation (in rates() because is pH dependent)
  conc_fresh['NH3'] <- (1/(1 + 10^(- log_ka['NH3'] - pH))) * conc_fresh[['TAN']] # Free NH3 concentration (g N/kg)
  conc_fresh['NH4'] <- conc_fresh[['TAN']] - conc_fresh[['NH3']]                 # (g N/kg)
  H2S_frac <- (1 - (1/(1 + 10^(- log_ka[['H2S']] - pH))))                        # H2S fraction of total sulfide
  pH_inhib <- (1 + 2*10^(0.5*(pH_lwr - pH_upr)))/(1 + 10^(pH - pH_upr) + 10^(pH_lwr - pH)) # pH inhibition factor
  NH3_inhib <- ifelse(conc_fresh[['NH3']] <= ki_NH3_min, 1, exp(-2.77259*((conc_fresh[['NH3']] - ki_NH3_min)/(ki_NH3_max - ki_NH3_min))^2))
  NH4_inhib <- ifelse(conc_fresh[['NH4']] <= ki_NH4_min, 1, exp(-2.77259*((conc_fresh[['NH4']] - ki_NH4_min)/(ki_NH4_max - ki_NH4_min))^2))

  # H2S inhibition
  H2S_inhib_meth <- 1 - (H2S_frac*(sulfide/slurry_mass))/ki_H2S_meth # H2S _inhib factor for methanogens
  H2S_inhib_sr <- 1 - (H2S_frac*(sulfide/slurry_mass))/ki_H2S_sr # H2S _inhib factor for sr
  if (H2S_frac*(sulfide/slurry_mass)>ki_H2S_meth) H2S_inhib_meth <-0 
  if (H2S_frac*(sulfide/slurry_mass)>ki_H2S_sr) H2S_inhib_sr <-0 
  
  # rut has same elements as qhat--set here, calculated below
  rut <- NA * qhat
  
  # VFA consumption rate by sulfate reducers (g/d)
  rut[i_sr] <- ((qhat[i_sr] * VFA / slurry_mass * xa[i_sr] / slurry_mass) / (ks[i_sr] + VFA / slurry_mass)) * slurry_mass *
    (sulfate / slurry_mass) / (ks_SO4 + sulfate / slurry_mass) *
    NH3_inhib[i_sr] * NH4_inhib[i_sr] * H2S_inhib_sr * pH_inhib[i_sr] 

  # Substrate utilization rate by methanogen groups in g/d affected by inhibition terms
  rut[i_meth] <- ((qhat[i_meth] * VFA / slurry_mass * xa[i_meth] / slurry_mass) / (ks[i_meth] + VFA / slurry_mass) * 
                  slurry_mass) *
                 NH3_inhib[i_meth] * NH4_inhib[i_meth] * H2S_inhib_meth * pH_inhib[i_meth] 

  # H2S emission
  H2SEmissionRate <- (kl['H2S'] * area * H2S_frac * sulfide)

  # Sp consumption by aerobic respiration (g/d)
  respiration <- (kl['oxygen'] * area * ((kH_oxygen * 0.208) - 0))

  # Some checks for safety
  if (any(rut < 0)) stop('In rates() function rut < 0 or otherwise strange. Check qhat parameters (92gg7)')

  # Derivatives, all in g/d except slurry_mass = kg/d
  # NTS: Some of these repeated calculations could be moved up
  derivatives <- c(
    xa = yield * rut + xa_fresh * slurry_prod_rate - decay_rate * xa, # expands to multiple elements with element for each mic group
    slurry_mass = slurry_prod_rate,
    Sp = slurry_prod_rate * conc_fresh[['Sp']] - alpha * Sp + sum(decay_rate * xa) - respiration,
    VFA = alpha * Sp - sum(rut) + slurry_prod_rate * conc_fresh[['VFA']],
    sulfate = slurry_prod_rate * conc_fresh_SO4 - rut[i_sr] * COD_conv[['S']],
    sulfide = slurry_prod_rate * conc_fresh[['S2']] + rut[i_sr] * COD_conv[['S']] - H2SEmissionRate,
    CH4_emis_cum = sum(rut[i_meth]) * COD_conv[['CH4']],
    CO2_emis_cum =  sum(rut[i_meth]) * COD_conv[['CO2_anaer']] + rut[i_sr] * COD_conv[['CO2_sr']] + respiration * COD_conv[['CO2_aer']],
    COD_conv_cum = sum(rut[i_meth]) + respiration + rut[i_sr],
    COD_conv_cum_meth = sum(rut[i_meth]),
    COD_conv_cum_respir = respiration,
    COD_conv_cum_sr = rut[i_sr]
  )

  return(list(derivatives, c(NH4 = conc_fresh[['NH4']], NH3 = conc_fresh[['NH3']], 
                             H2S_inhib_meth = H2S_inhib_meth, H2S_inhib_sr = H2S_inhib_sr, 
                             pH_inhib = pH_inhib, qhat = qhat, NH3_inhib = NH3_inhib, 
                             NH4_inhib = NH4_inhib, rut = rut, respiration = respiration)))

}

rates <- function(t, 
		  y, 
		  parms, 
		  temp_C_fun = temp_C_fun, 
		  pH_fun = pH_fun, 
                  SO4_inhibition_fun = SO4_inhibition_fun) {

    y[y < 1E-10] <- 1E-10

    # Remove slurry mass from parms to not overwrite y['slurry_mass']
    parms$slurry_mass <- NULL
    
    # all elements in parms exported to current environment
    list2env(parms, envir = environment())

    pH <- pH_fun(t + t_run)
   
    temp_C <- temp_C_fun(t + t_run)
    temp_K <- temp_C + 273.15

    # Extract state variable values from y argument
    xa <- y[1:n_mic]

    # Move elements of y into rates environment
    list2env(as.list(y[-c(1:n_mic)]), envir = environment())

    # Hydrolysis rate . . . 
    alpha <- 0.1    

    # Microbial substrate utilization rate (vectorized calculation)
    qhat <- scale[['qhat_opt']] * CTM_cpp(temp_K, T_opt, T_min, T_max, qhat_opt)
    names(qhat) <- names(qhat_opt)
    
    # Decay of all microorganisms follow same kinetics with faster decay at higher temp up until 313, at which constant decay of 0.02
    
    decay_rate <- CTM_cpp(temp_K, 313, 273, 325, decay_rate)
    if(temp_K > 313) decay_rate <- 0.02
    
    
    rates_inhib <- combined_cpp(g_NH4, pH_inhib_overrule,
      pH, pH_floor,
      TAN, VFA, sulfide, slurry_mass,
      pH_LL, pH_UL,
      ki_NH3_min, ki_NH3_max,
      ki_NH4_min, ki_NH4_max,
      ki_HAC, ki_H2S_slope, 
      ki_H2S_int, ki_H2S_min,
      IC50_low, 
      temp_K, temp_C, temp_standard, 
      area, floor_area,
      Cfat, CPs, CPf, RFd, 
      starch, VSd, resp, kl,
      qhat, i_meth, i_sr, xa,
      ks_coefficient, scale_ks = scale[['ks_coefficient']], ks_SO4, sulfate,
      urea, alpha_urea = alpha[['urea']], km_urea)
    
    list2env(rates_inhib, envir = environment())
  
    # CO2 production from fermentation + methanogenesis + sulfate reduction + aerobic respiration at slurry surface
    # also calcualtes growth rate of xa_bac and xa_aer, mineralization rates and COD production from fermentation
    
    ### MOVE all the stuff below until derivatives to CPP?
    ferm <- stoich(alpha, y, conc_fresh, sub_resp, respiration,
                   carb, pro, lip, carb_resp, pro_resp, lip_resp, ace, hyd,
                   ace_sr, hyd_sr)
   
    ferm$COD_conv_meth_CO2 <- COD_conv_meth_CO2
    ferm$COD_conv_sr_CO2 <- COD_conv_sr_CO2
    ferm$TAN_min_ferm <- TAN_min_ferm
    ferm$TAN_min_resp <- TAN_min_resp * respiration
    ferm$xa_aer_rate <- xa_aer_rate * respiration
    ferm$xa_bac_rate <- xa_bac_rate
    ferm$VFA_H2 <- VFA_H2_ferm
    ferm$CO2_resp <- CO2_resp * respiration
    ferm$ferm[['CO2']] <- CO2_ferm

    
    CO2_ferm_meth_sr <- ferm$ferm[['CO2']] * 44.01 + sum(rut[i_meth+1])/ferm$COD_conv_meth_CO2[[1]] + sum(rutsr)/ferm$COD_conv_sr_CO2[[1]]
    
    # Derivatives, all in gCOD/d except slurry_mass = kg/d, N and S are gN or gS, VSd_A and VSnd_A are g/d

    derivatives <- derivatives_cpp(
      yield, rut, slurry_prod_rate, decay_rate, 
      scale, xa_fresh, conc_fresh, 
      ferm, alpha, COD_conv, 
      rain, evap,  area,  respiration,  sub_resp,
       rut_urea,  NH3_emis_rate_pit,  NH3_emis_rate_floor, 
       N2O_emis_rate, rutsr,  H2S_emis_rate, 
       R,  temp_K, i_meth, i_sr,  CO2_ferm_meth_sr, 
      xa,  slurry_mass,  xa_bac,  xa_aer,  xa_dead,
       RFd,  iNDF,  ash,  VSd,  starch,  CPs,  CPf,
       Cfat,  VFA,  urea,  TAN,  sulfate,  sulfide, qhat
    )

    names(derivatives) <- names(y)

    return(list(derivatives, c(COD_load_rate = derivatives[['COD_load_cum']], C_load_rate = derivatives[['C_load_cum']], N_load_rate = derivatives[['N_load_cum']],
                               CH4_emis_rate = derivatives[['CH4_emis_cum']], CO2_emis_rate = derivatives[['CO2_emis_cum']], 
                               H2S_emis_rate = H2S_emis_rate, NH3_emis_rate_pit = NH3_emis_rate_pit, NH3_emis_rate_floor = NH3_emis_rate_floor,
                               N2O_emis_rate = N2O_emis_rate,
                               qhat = qhat, alpha = alpha, H2S_inhib = H2S_inhib, NH3_inhib = NH3_inhib, NH4_inhib = NH4_inhib, HAC_inhib = HAC_inhib, cum_inhib = cum_inhib, 
                               conc_fresh = conc_fresh, xa_init = xa_init, 
                               xa_fresh = xa_fresh * scale[['xa_fresh']], area = area, slurry_prod_rate = slurry_prod_rate,
                               respiration = respiration, rain = rain, evap = evap, CO2_ferm_meth_sr = CO2_ferm_meth_sr)))
   
}

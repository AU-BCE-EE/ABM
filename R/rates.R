rates <- function(t, y, parms, temp_C_fun = temp_C_fun, pH_fun = pH_fun, 
                  SO4_inhibition_fun = SO4_inhibition_fun, 
                  conc_fresh_fun = conc_fresh_fun, xa_fresh_fun = xa_fresh_fun) {

    y[y < 1E-10] <- 1E-10

    # need to remove slurry mass from parms to not overwrite y['slurry_mass']
    parms$slurry_mass <- NULL
    
    # all elements in parms exported to current environment
    list2env(parms, envir = environment())

    # correct slurry production rate in periods with grazing
    if(!is.null(graze_int) & any(graze_int != 0)){
       slurry_prod_rate <- graze_fun(t,  t_run, days, slurry_prod_rate, graze_int, graze_hours = graze[['hours_day']])
    }

    # pH, numeric, variable, or from H2SO4
    if (is.numeric(pH) | is.data.frame(pH)) {
      pH <- pH_fun(t + t_run)
    } else if (pH == 'calc') {
      pH <- H2SO4_titrat(conc_SO4 = conc_fresh[['sulfate']], class_anim = "pig")$pH
    } else {
      stop('pH problem (xi342)')
    }
    # Remove name 'pH' that sometimes comes along
    pH <- as.numeric(pH)
    
    # For time-variable fresh concentrations, need to use function to get fresh concentrations at particular time
    if(is.data.frame(conc_fresh)){
      conc_fresh <- list()
      conc_fresh$sulfide <- conc_fresh_fun$conc_fresh_fun_sulfide(t + t_run)
      conc_fresh$urea <- conc_fresh_fun$conc_fresh_fun_urea(t + t_run)
      conc_fresh$sulfate <- conc_fresh_fun$conc_fresh_fun_sulfate(t + t_run)
      conc_fresh$TAN <- conc_fresh_fun$conc_fresh_fun_TAN(t + t_run)
      conc_fresh$starch <- conc_fresh_fun$conc_fresh_fun_starch(t + t_run)
      conc_fresh$VFA <- conc_fresh_fun$conc_fresh_fun_VFA(t + t_run)
      conc_fresh$xa_aer <- conc_fresh_fun$conc_fresh_fun_xa_aer(t + t_run)
      conc_fresh$xa_bac <- conc_fresh_fun$conc_fresh_fun_xa_bac(t + t_run)
      conc_fresh$xa_dead <- conc_fresh_fun$conc_fresh_fun_xa_dead(t + t_run)
      conc_fresh$Cfat <- conc_fresh_fun$conc_fresh_fun_Cfat(t + t_run)
      conc_fresh$CPs <- conc_fresh_fun$conc_fresh_fun_CPs(t + t_run)
      conc_fresh$CPf <- conc_fresh_fun$conc_fresh_fun_CPf(t + t_run)
      conc_fresh$RFd <- conc_fresh_fun$conc_fresh_fun_RFd(t + t_run)
      conc_fresh$iNDF <- conc_fresh_fun$conc_fresh_fun_iNDF(t + t_run)
      conc_fresh$VSd <- conc_fresh_fun$conc_fresh_fun_VSd(t + t_run)
      conc_fresh$VSd_A <- conc_fresh_fun$conc_fresh_fun_VSd_A(t + t_run)
      conc_fresh$VSnd_A <- conc_fresh_fun$conc_fresh_fun_VSnd_A(t + t_run)
      conc_fresh$ash <- conc_fresh_fun$conc_fresh_fun_ash(t + t_run)
    }
    
    #INTERPOLATE CPP. Check this
    if(is.data.frame(xa_fresh)) {
      xa_fresh <- sapply(seq_along(grps), function(i) {
        conc <- xa_fresh_fun[[i]](t + t_run)
        return(conc)
      })
      names(xa_fresh) <- grps
    }

    #temp functions
    # CPP interpolation func?
    temp_C <- temp_C_fun(t + t_run)
    temp_K <- temp_C + 273.15

    # Extract state variable values from y argument
    xa <- y[1:n_mic]

    # Move elements of y into rates environment
    list2env(as.list(y[-c(1:n_mic)]), envir = environment())

    # Hydrolysis rate with Arrhenius function cpp. 
    alpha <-  Arrh_func_cpp(A, E, R, temp_K, scale_alpha_opt = scale[['alpha_opt']], alpha_opt_scale_type, alpha_opt_scale_CP)
    names(alpha) <- names(A)
    
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
    
    CO2_ferm_meth_sr <- ferm$ferm[['CO2']] * 44.01 + sum(rut[i_meth+1])/ferm$COD_conv_meth_CO2[[1]] + sum(rutsr)/ferm$COD_conv_sr_CO2[[1]]
    CO2_ferm <- ferm$ferm[['CO2']] * 44.01 
    
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
                               respiration = respiration, rain = rain, evap = evap)))
   
}

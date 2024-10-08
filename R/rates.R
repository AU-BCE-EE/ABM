rates <- function(t, y, parms, temp_C_fun = temp_C_fun, pH_fun = pH_fun, 
                  SO4_inhibition_fun = SO4_inhibition_fun, 
                  conc_fresh_fun = conc_fresh_fun, xa_fresh_fun = xa_fresh_fun) {

    y[y < 1E-10] <- 1E-10
     
    # need to remove slurry mass from parms to not overwrite y['slurry_mass']
    parms$slurry_mass <- NULL
     
    # Put all parameters in parms elements directly in rates environment
    for (pp in names(parms)) {
      assign(pp, parms[[pp]])
    }
  
    # correct slurry production rate in periods with grazing
    suppressWarnings({
    if(!is.null(graze_int) & any(graze_int) != 0) {
       slurry_prod_rate <- graze_fun(t,  t_run, days, slurry_prod_rate, graze_int, graze_hours = graze[['hours_day']])
    }
    })
       
    # pH, numeric, variable, or from H2SO4
    if (is.numeric(pH) | is.data.frame(pH)) {
      pH <- pH_fun(t + t_run)
    } else if (pH == 'calc') {
      pH <- H2SO4_titrat(conc_SO4 = conc_fresh[['sulfate']], class_anim = "pig")$pH
    } else {
      stop('pH problem (xi342)')
    }
    
    # calculate time of a batch
    if (!is.na(wash_int)){
      batches <- c(floor((t + t_run)/(wash_int + rest_d)))
      t_batch <- (t + t_run) - batches * (wash_int + rest_d)
      if (t_batch > wash_int) t_batch <- 0
    } else {
      t_batch <- 0
    }
    
    # if urea fresh increase during a batch 
    if (!is.na(slopes[['urea']]) & !is.data.frame(conc_fresh)) {
      start_urea <- conc_fresh[['urea']] - slopes[['urea']] * wash_int/2
      conc_fresh[['urea']] <- slopes[['urea']] * t_batch + start_urea
    }
    
    # if slurry production increases during a batch
    slurry_prod_rate_default <- slurry_prod_rate
    
    if (!is.na(slopes[['slurry_prod_rate']]) && slurry_prod_rate_default != 0 && !is.na(wash_int)) {
      start_slurry_prod_rate <- slurry_prod_rate_default - slopes[['slurry_prod_rate']] * wash_int/2
      slurry_prod_rate <- slopes[['slurry_prod_rate']] * t_batch + start_slurry_prod_rate
      if (t_batch > wash_int) slurry_prod_rate <- 0
    }
    
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
      conc_fresh$CP <- conc_fresh_fun$conc_fresh_fun_CP(t + t_run)
      conc_fresh$RFd <- conc_fresh_fun$conc_fresh_fun_RFd(t + t_run)
      conc_fresh$iNDF <- conc_fresh_fun$conc_fresh_fun_iNDF(t + t_run)
      conc_fresh$VSd <- conc_fresh_fun$conc_fresh_fun_VSd(t + t_run)
      conc_fresh$VSd_A <- conc_fresh_fun$conc_fresh_fun_VSd_A(t + t_run)
      conc_fresh$VSnd_A <- conc_fresh_fun$conc_fresh_fun_VSnd_A(t + t_run)
      conc_fresh$ash <- conc_fresh_fun$conc_fresh_fun_ash(t + t_run)
    }
    
    xa_fresh <- if (is.data.frame(xa_fresh)) {
      sapply(seq_along(grps), function(i) {
        conc <- xa_fresh_fun[[i]](t + t_run)
        return(conc)
      })
    } else{
      xa_fresh <- xa_fresh
    }
    
    names(xa_fresh) <- grps
    
    # Hard-wired temp settings settings
    temp_standard <- 298
    temp_zero <- 273
    
    #temp functions
    temp_C <- temp_C_fun(t + t_run)
    temp_K <- temp_C + 273.15

    
    # Find methanogens and sulfate reducers
    i_meth <- grepl('^[mp]', names(qhat_opt))
    i_sr <- grepl('^sr', names(qhat_opt))
    n_mic <- length(qhat_opt)
    
    # Extract state variable values from y argument
    xa <- y[1:n_mic]
    y <- as.list(y[-c(1:n_mic)])

    # Move elements of y into rates environment
    for (pp in names(y)) {
      assign(pp, y[[pp]])
    }
    
    # Hard-wired equilibrium constants
    log_ka <- c(NH3 = - 0.09046 - 2729.31/temp_K, 
                H2S = - 3448.7/temp_K + 47.479 - 7.5227* log(temp_K),
                HAC = -4.8288 + 21.42/temp_K)
    
    kH_oxygen <- 0.0013 * exp(1700 * ((1 / temp_K) - (1 / temp_standard))) * 32 * 1000
   
    # Hard-wire NH4+ activity coefficient
    g_NH4 <- 0.7

    # Hydrolysis rate with Arrhenius function or CTM. 
    alpha <-  Arrh_func(A, E, R, temp_K)

    if(conc_fresh$VSd <= 10E-9){
      alpha_opt_scale_type <- scale_alpha_opt[['notVSd']]
      alpha_opt_scale_CP <- scale_alpha_opt[['CP']]
    } else if(conc_fresh$VSd > 10E-9){
      alpha_opt_scale_type  <- scale_alpha_opt[['VSd']]
      alpha_opt_scale_CP <- 1
    }
    
    alpha[names(alpha) != 'urea'] <- scale[['alpha_opt']] * alpha_opt_scale_type * alpha[names(alpha) != 'urea']
    alpha['CP'] <- alpha_opt_scale_CP * alpha['CP']
    
    # Microbial substrate utilization rate (vectorized calculation)
    qhat <- scale[['qhat_opt']] * CTM_cpp(temp_K, T_opt, T_min, T_max, qhat_opt)
    names(qhat) <- names(qhat_opt)
    
    # Ks temperature dependence
    ks <- ks_coefficient * (0.8157 * exp(-0.063 * temp_C)) 
    
    # NTS: Move this all out to a speciation function???
    # Chemical speciation (in rates() because is pH dependent)
    pH_surf <- pH + 1 # rough approximation from Rotz et al. IFSM 2012., "markfoged" thesis, "Petersen et al. 2014", "bilds?e et al. not published", "Elzing & Aarnik 1998", and own measurements..
    pH_surf_floor <- 8 # pH of manure on floor is not affected by acidification. Kept at 6.8
    
    HAC_frac <- 1-(1/(1 + 10^(-log_ka[['HAC']] - pH)))
    NH3_frac <- ((1/(1 + 10^(- log_ka[['NH3']] + log10(g_NH4) - pH)))) 
    NH3_frac_surf <- ((1/(1 + 10^(- log_ka[['NH3']] + log10(g_NH4) - pH_surf))))
    NH3_frac_floor <- ((1/(1 + 10^(- log_ka[['NH3']] + log10(g_NH4) - pH_surf_floor))))
    H2S_frac <- 1 - (1/(1 + 10^(- log_ka[['H2S']] - pH))) # H2S fraction of total sulfide

    # NTS: or just add H2CO3* here?
    # NTS: need TIC production too
    
    pH_inhib <- 1
    NH3_inhib <- 1
    NH4_inhib <- 1
    HAC_inhib <- 1
    H2S_inhib <- 1
    
    # if pH_inhibiion should be used, NH4 and NH3 inhibition is ignored and pH inhibition is used instead. 
    if(pH_inhib_overrule){
            pH_inhib <- (1 + 2*10^(0.5* (pH_LL - pH_UL)))/(1+ 10^(pH - pH_UL) + 10^(pH_LL - pH))
    } else{
      #NH3 and NH4 inhibition      
      NH3_inhib <- ifelse(NH3_frac * TAN/slurry_mass <= ki_NH3_min, 1, exp(-2.77259 * ((NH3_frac * (TAN/(slurry_mass)) - ki_NH3_min)/(ki_NH3_max - ki_NH3_min))^2))
      NH4_inhib <- ifelse((1 - NH3_frac) * (TAN/slurry_mass) <= ki_NH4_min, 1, exp(-2.77259*(((1 - NH3_frac) * (TAN/(slurry_mass)) - ki_NH4_min)/(ki_NH4_max - ki_NH4_min))^2))
      # HAC inhibition
      HAC_inhib <- ifelse(HAC_frac * VFA/slurry_mass >= 0.05, (2-ki_HAC/(ki_HAC+0.05)) * ki_HAC/(ki_HAC + (HAC_frac * VFA/slurry_mass)), 1) 
      # H2S inhibition
      H2S_inhib <- NA * qhat
            
      if(pH >= 6.8){
        IC50 <- ki_H2S_slope * pH + ki_H2S_int
      } else {
        IC50 <- IC50_low
      } 
      
      a <- -0.5/(IC50 - (H2S_frac * ki_H2S_min))
      x <- H2S_frac * sulfide/(slurry_mass)
      b <- 1 -(-0.5/(IC50 - (H2S_frac * ki_H2S_min)) * H2S_frac * ki_H2S_min/(slurry_mass))
      
      H2S_inhib <- a * x + b
      H2S_inhib[H2S_inhib < 0] <- 0
      H2S_inhib[H2S_inhib > 1 ] <- 1
    
    }
    
    cum_inhib <- HAC_inhib * NH3_inhib * NH4_inhib * H2S_inhib * pH_inhib
      
    # Henrys constant temp dependency
    H.NH3 <- 1431 * 1.053^(293 - temp_K)
    
    # Reduction from cover 
    kl[['NH3']] <- kl[['NH3']] * EF_NH3 
  
    # NH3 emission g(N) pr day
    NH3_emis_rate_floor <- kl[['NH3_floor']] * floor_area * ((NH3_frac_floor * TAN)/(slurry_mass)) * 1000 / H.NH3 # multiply by 1000 to get from g/kg to g/m3
    NH3_emis_rate_pit <- kl[['NH3']] * area * ((NH3_frac_surf * TAN)/(slurry_mass)) * 1000 / H.NH3 # multiply by 1000 to get from g/kg to g/m3
    
    # N2O emission g(N) pr day
    N2O_emis_rate <- area * EF_N2O # 
    
    # H2S emission g(S) pr day
    H2S_emis_rate <- kl[['H2S']] * area * ((H2S_frac * sulfide)/(slurry_mass)) * 1000 # multiply by 1000 to get from g/kg to g/m3

    # Respiration gCOD/d, second-order reaction where kl applies for substrate concentration of 100 g COD / kg slurry
    respiration <- 0
    sub_resp <- Cfat + CP + RFd + starch + VSd 
    
    # Solver is slowed significantly when the slurry_mass is low relative to surface area is low, therefore cut off at 1 mm slurry height
    if(resp & slurry_mass/area >= 1){
      kl.oxygen <- exp(0.6158816 + 0.09205127 * temp_C) # from own lab experiments (Dalby et al. 2024..unpublished) 
      respiration <- kl.oxygen * area * ((kH_oxygen * 0.208) - 0) * (sub_resp / slurry_mass) / 100
    }
    
    # VFA uptake rates
    rut <- NA * qhat
    
    # VFA consumption rate by sulfate reducers (g/d) affected by inhibition terms
    rut[i_sr] <- ((qhat[i_sr] * VFA / (slurry_mass) * xa[i_sr] / (slurry_mass) / (scale[['ks_coefficient']] * ks[i_sr] + VFA / (slurry_mass))) * (slurry_mass) *
                    (sulfate / (slurry_mass)) / (ks_SO4 + sulfate / (slurry_mass))) * cum_inhib[i_sr]

    # VFA consumption rate by methanogen groups (g/d) affected by inhibition terms
    rut[i_meth] <- ((qhat[i_meth] * VFA / (slurry_mass) * xa[i_meth] / (slurry_mass)) / (scale[['ks_coefficient']] * ks[i_meth] + VFA / (slurry_mass)) *
                      (slurry_mass)) * cum_inhib[i_meth]
    
    # urea hydrolysis by Michaelis Menten
    rut_urea <- (alpha[['urea']] * urea) / (km_urea + urea/slurry_mass)
    
    # Some checks for safety
    if (any(rut < 0)) stop('In rates() function rut < 0 or otherwise strange. Check qhat parameters (92gg7)')

    # If there are no SRs...
    rutsr <- rut[i_sr]
    if (length(rutsr) == 0) rutsr <- 0

    # CO2 production from fermentation + methanogenesis + sulfate reduction + aerobic respiration at slurry surface
    # also calcualtes growth rate of xa_bac and xa_aer, mineralization rates and COD production from fermentation
    
    ferm <- stoich(alpha, y, conc_fresh, sub_resp, respiration)
    CO2_ferm_meth_sr <- ferm$ferm['CO2'] * 44.01 + sum(rut[i_meth])/ferm$COD_conv_meth_CO2 + sum(rutsr)/ferm$COD_conv_sr_CO2
    CO2_ferm <- ferm$ferm['CO2'] * 44.01 
    
    # Derivatives, all in gCOD/d except slurry_mass = kg/d, N and S are gN or gS, VSd_A and VSnd_A are g/d
    # NTS: Some of these repeated calculations could be moved up
    # need to implement xa_bac to keep mass balance here: 

    derivatives <- c(
       xa = scale[['yield']] * yield * rut + scale[['xa_fresh']] * xa_fresh * slurry_prod_rate - decay_rate * xa, # expands to multiple elements with element for each mic group
       slurry_mass = slurry_prod_rate + (rain - evap) * area,
       xa_aer = slurry_prod_rate * conc_fresh[['xa_aer']] + ferm[['xa_aer_rate']] - decay_rate_xa * xa_aer,
       xa_bac = slurry_prod_rate * conc_fresh[['xa_bac']] + ferm[["xa_bac_rate"]] - decay_rate_xa * xa_bac, # xa_bac_rate is the growth pr day (calculated in stoich function, therefore not necessary to multiply with a yield coeff)
       xa_dead = slurry_prod_rate * conc_fresh[['xa_dead']] - alpha[['xa_dead']] * xa_dead + sum(decay_rate * xa) + decay_rate_xa * (xa_bac + xa_aer) * (1 - COD_conv[['frac_CP_xa']]),
       RFd = slurry_prod_rate * conc_fresh[['RFd']] - alpha[['RFd']] * RFd - respiration * RFd/sub_resp,
       iNDF = slurry_prod_rate * conc_fresh[['iNDF']],
       ash = slurry_prod_rate * conc_fresh[['ash']],
       VSd = slurry_prod_rate * conc_fresh[['VSd']] - alpha[['VSd']] * VSd - respiration * VSd/sub_resp,
       starch = slurry_prod_rate * conc_fresh[['starch']] - alpha[['starch']] * starch - respiration * starch/sub_resp,
       CP = slurry_prod_rate * conc_fresh[['CP']] - alpha[['CP']] * CP - respiration * CP/sub_resp + decay_rate_xa * (xa_bac + xa_aer) * COD_conv[['frac_CP_xa']],
       Cfat = slurry_prod_rate * conc_fresh[['Cfat']] - alpha[['Cfat']] * Cfat - respiration * Cfat/sub_resp,
       VFA = alpha[['xa_dead']] * xa_dead + ferm[["VFA_H2"]] - sum(rut) + slurry_prod_rate * conc_fresh[['VFA']],
       urea = slurry_prod_rate * conc_fresh[['urea']] - rut_urea,
       TAN = slurry_prod_rate * conc_fresh[['TAN']] + rut_urea + ferm[['TAN_min_ferm']] + ferm[['TAN_min_resp']] - NH3_emis_rate_pit - NH3_emis_rate_floor - N2O_emis_rate,
       sulfate = slurry_prod_rate * conc_fresh[['sulfate']] - sum(rutsr) / COD_conv[['S']],
       sulfide = slurry_prod_rate * conc_fresh[['sulfide']] + sum(rutsr) / COD_conv[['S']] - H2S_emis_rate,
       VSd_A = - VSd_A * ((exp(lnA[['VSd_A']] - E_CH4[['VSd_A']] / (R * temp_K))) * 24 / 1000 * VS_CH4) + slurry_prod_rate * conc_fresh[['VSd_A']],
       VSnd_A = slurry_prod_rate * conc_fresh[['VSnd_A']],
       CH4_A_emis_cum = VSd_A * (exp(lnA[['VSd_A']] - E_CH4[['VSd_A']] / (R * temp_K))) * 24 / 1000,
       NH3_emis_cum = NH3_emis_rate_pit + NH3_emis_rate_floor,
       N2O_emis_cum = N2O_emis_rate,
       CH4_emis_cum = sum(rut[i_meth]) / COD_conv[['CH4']],
       CO2_emis_cum = CO2_ferm_meth_sr + ferm[['CO2_resp']] + rut_urea / COD_conv[['CO2_ureo']],
       COD_conv_cum = sum(rut[i_meth]) + respiration + sum(rutsr),
       COD_conv_cum_meth = sum(rut[i_meth]),
       COD_conv_cum_respir = respiration,
       COD_conv_cum_sr = rutsr
     )

    return(list(derivatives, c(H2S_emis_rate = H2S_emis_rate, NH3_emis_rate_pit = NH3_emis_rate_pit,
                               NH3_emis_rate_floor = NH3_emis_rate_floor,
                               qhat = qhat, alpha = alpha, CO2_ferm = CO2_ferm, CO2_ferm_meth_sr = CO2_ferm_meth_sr, CO2_resp = ferm[['CO2_resp']],
                               H2S_inhib = H2S_inhib, NH3_inhib = NH3_inhib, NH4_inhib = NH4_inhib,
                               TAN_min_resp = ferm[['TAN_min_resp']], TAN_min_ferm = ferm[['TAN_min_ferm']], HAC_inhib = HAC_inhib, cum_inhib = cum_inhib, 
                               rut = rut, rut_urea = rut_urea, t_run = t_run, t_batch = t_batch, conc_fresh = conc_fresh, xa_init = xa_init, 
                               xa_fresh = xa_fresh * scale[['xa_fresh']], area = area, t_batch = t_batch, slurry_prod_rate = slurry_prod_rate,
                               respiration = respiration, rain = rain, evap = evap)))
}

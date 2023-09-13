rates <- function(t, y, parms, temp_C_fun = temp_C_fun, pH_fun = pH_fun, 
                  SO4_inhibition_fun = SO4_inhibition_fun, 
                  conc_fresh_fun = conc_fresh_fun, xa_fresh_fun = xa_fresh_fun) {
  
     y[y < 1E-10] <- 1E-10
     
     # need to remove slurry mass from parms to not overwrite y['slurry_mass']
     parms$slurry_mass <- NULL
     
     with(as.list(parms), {
    
    # correct slurry production rate in periods with grazing
    if(!is.null(graze_int) & graze_int != 0) {
       slurry_prod_rate <- graze_fun(t,  t_run, days, slurry_prod_rate, graze_int, graze_hours)
    }
       
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
    
    if (!is.na(slopes[['slurry_prod_rate']]) && slurry_prod_rate_default != 0) {
      start_slurry_prod_rate <- slurry_prod_rate_default - slopes[['slurry_prod_rate']] * wash_int/2
      slurry_prod_rate <- slopes[['slurry_prod_rate']] * t_batch + start_slurry_prod_rate
      if (t_batch > wash_int) slurry_prod_rate <- 0
    }

    # when variable fresh concentration is used
    #if(is.data.frame(conc_fresh) | is.data.frame(xa_fresh)){
    #  conc_fresh <- variable_conc(conc_fresh, xa_fresh, t, t_run)$conc_fresh
    #  xa_fresh <- variable_conc(conc_fresh, xa_fresh, t, t_run)$xa_fresh
    #}
    
    if(is.data.frame(conc_fresh)){
      conc_fresh <- list()
      conc_fresh$sulfide <- conc_fresh_fun$conc_fresh_fun_sulfide(t + t_run)
      conc_fresh$urea <- conc_fresh_fun$conc_fresh_fun_sulfide(t + t_run)
      conc_fresh$sulfate <- conc_fresh_fun$conc_fresh_fun_sulfate(t + t_run)
      conc_fresh$TAN <- conc_fresh_fun$conc_fresh_fun_TAN(t + t_run)
      conc_fresh$starch <- conc_fresh_fun$conc_fresh_fun_starch(t + t_run)
      conc_fresh$VFA <- conc_fresh_fun$conc_fresh_fun_VFA(t + t_run)
      conc_fresh$xa_dead <- conc_fresh_fun$conc_fresh_fun_xa_dead(t + t_run)
      conc_fresh$CF <- conc_fresh_fun$conc_fresh_fun_CF(t + t_run)
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
    list2env(y, envir = .GlobalEnv)
    
    # Hard-wired equilibrium constants
    log_ka <- c(NH3 = - 0.09046 - 2729.31/temp_K, 
                H2S = - 3448.7/temp_K + 47.479 - 7.5227* log(temp_K),
                HAC = -4.8288 + 21.42/temp_K)
    
    kH_oxygen <- 0.0013 * exp(1700 * ((1 / temp_K) - (1 / temp_standard))) * 32 * 1000
   
    # Hard-wire NH4+ activity coefficient
    g_NH4 <- 0.7
    
    # Hydrolysis rate with Arrhenius function or CTM. 
    alpha <- scale['alpha_opt'] * CTM_cpp(temp_K, alpha_T_opt, alpha_T_min, alpha_T_max, alpha_opt)
    names(alpha) <- names(alpha_opt)
    alpha_arrh <- Arrh_func(A, E, R, temp_K)
    alpha <- c(alpha, alpha_arrh)

    # Microbial substrate utilization rate (vectorized calculation)
    qhat <- scale['qhat_opt'] * CTM_cpp(temp_K, T_opt, T_min, T_max, qhat_opt)
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
    
    # HAC inhibition
    HAC_inhib <- ifelse(HAC_frac * VFA/slurry_mass >= 0.05, 1.16*0.31/(0.31 + (HAC_frac * VFA/slurry_mass)), 1) 

    # NH3, NH4 inhibition
    NH3_inhib <- ifelse(NH3_frac * TAN/slurry_mass <= ki_NH3_min, 1, exp(-2.77259 * ((NH3_frac * (TAN/(slurry_mass)) - ki_NH3_min)/(ki_NH3_max - ki_NH3_min))^2))
    NH4_inhib <- ifelse((1 - NH3_frac) * (TAN/slurry_mass) <= ki_NH4_min, 1, exp(-2.77259*(((1 - NH3_frac) * (TAN/(slurry_mass)) - ki_NH4_min)/(ki_NH4_max - ki_NH4_min))^2))
  
    # H2S inhibition
    H2S_inhib <- NA * qhat
    
        if(pH >= 6.8){
          IC50 <- ki_H2S_slope * pH + ki_H2S_int
        } else {
          IC50 <- ki_H2S_slope * 6.8 + ki_H2S_int
        } 
      
        a <- -0.5/(IC50 - (H2S_frac * ki_H2S_min))
        x <- H2S_frac * sulfide/(slurry_mass)
        b <- 1 -(-0.5/(IC50 - (H2S_frac * ki_H2S_min)) * H2S_frac * ki_H2S_min/(slurry_mass))
      
        H2S_inhib <- a * x + b
        H2S_inhib[H2S_inhib < 0] <- 0
        H2S_inhib[H2S_inhib > 1 ] <- 1
      
      cum_inhib <- HAC_inhib * NH3_inhib * NH4_inhib * H2S_inhib
      
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
    
    kl.oxygen <- (0.1699211 * temp_C) # from own lab experiments (Dalby et al. 2023..unpublished) 
    
    sub_respir <- CF + CP + RFd + starch + VSd 
    if (sub_respir <= 0) sub_respir <- 1E-20
    respiration <- kl.oxygen * area * ((kH_oxygen * 0.208) - 0) * (sub_respir / slurry_mass) / 100 
    
    # VFA uptake rates
    rut <- NA * qhat
    
    # VFA consumption rate by sulfate reducers (g/d) affected by inhibition terms
    rut[i_sr] <- ((qhat[i_sr] * VFA / (slurry_mass) * xa[i_sr] / (slurry_mass) / (scale['ks_coefficient'] * ks[i_sr] + VFA / (slurry_mass))) * (slurry_mass) *
                    (sulfate / (slurry_mass)) / (ks_SO4 + sulfate / (slurry_mass))) * cum_inhib[i_sr]

    # VFA consumption rate by methanogen groups (g/d) affected by inhibition terms
    rut[i_meth] <- ((qhat[i_meth] * VFA / (slurry_mass) * xa[i_meth] / (slurry_mass)) / (scale['ks_coefficient'] * ks[i_meth] + VFA / (slurry_mass)) *
                      (slurry_mass)) * cum_inhib[i_meth]
   
    # Some checks for safety
    if (any(rut < 0)) stop('In rates() function rut < 0 or otherwise strange. Check qhat parameters (92gg7)')

    # If there are no SRs...
    rutsr <- rut[i_sr]
    if (length(rutsr) == 0) rutsr <- 0

    # CO2 production from fermentation + methanogenesis + sulfate reduction
    ferm <- stoich(alpha, y, conc_fresh)
    CO2_ferm_meth_sr <- ferm$ferm['CO2'] * 44.01 + sum(rut[i_meth])/ferm$COD_conv_meth_CO2 + sum(rutsr)/ferm$COD_conv_sr_CO2
    CO2_ferm <- ferm$ferm['CO2'] * 44.01 
    
    # surface respiration conv_COD_aer assumes 10% cell yield, but aerobic biomass is not tracked. 
    # Therefore aer cell mass should be accounted for somewhere. Suggest to add to the xa_dead pool for COD balance 
    xa_aer <- respiration * 0.1
    
    # Derivatives, all in g/d except slurry_mass = kg/d
    # NTS: Some of these repeated calculations could be moved up
    derivatives <- c(
       xa = scale['yield'] * yield * rut + scale['xa_fresh'] * xa_fresh * slurry_prod_rate - decay_rate * xa, # expands to multiple elements with element for each mic group
       slurry_mass = slurry_prod_rate + (rain - evap) * area,
       xa_dead = slurry_prod_rate * conc_fresh[['xa_dead']] - alpha[['xa_dead']] * xa_dead + sum(decay_rate * xa) + xa_aer,
       RFd = slurry_prod_rate * conc_fresh[['RFd']] - alpha[['RFd']] * RFd - respiration * RFd/sub_respir,
       iNDF = slurry_prod_rate * conc_fresh[['iNDF']],
       ash = slurry_prod_rate * conc_fresh[['ash']],
       VSd = slurry_prod_rate * conc_fresh[['VSd']] - alpha[['VSd']] * VSd - respiration * VSd/sub_respir,
       starch = slurry_prod_rate * conc_fresh[['starch']] - alpha[['starch']] * starch - respiration * starch/sub_respir,
       CP = slurry_prod_rate * conc_fresh[['CP']] - alpha[['CP']] * CP - respiration * CP/sub_respir,
       CF = slurry_prod_rate * conc_fresh[['CF']] - alpha[['CF']] * CF - respiration * CF/sub_respir,
       VFA = alpha[['xa_dead']] * xa_dead + alpha[['starch']] * starch + alpha[['CP']] * CP + alpha[['CF']] * CF + alpha[['RFd']] * RFd + alpha[['VSd']] * VSd   - sum(rut) + slurry_prod_rate * conc_fresh[['VFA']],
       urea = slurry_prod_rate * conc_fresh[['urea']] - alpha[['urea']] * urea,
       TAN = slurry_prod_rate * conc_fresh[['TAN']] + alpha[['urea']] * urea + alpha[['CP']] * CP / COD_conv[['N_CP']] + respiration * CP/sub_respir / COD_conv[['N_CP']] - NH3_emis_rate_pit - NH3_emis_rate_floor - N2O_emis_rate,
       sulfate = slurry_prod_rate * conc_fresh[['sulfate']] - sum(rutsr) / COD_conv[['S']],
       sulfide = slurry_prod_rate * conc_fresh[['sulfide']] + sum(rutsr) / COD_conv[['S']] - H2S_emis_rate,
       VSd_A = - VSd_A * ((exp(lnA[['VSd_A']] - E_CH4[['VSd_A']] / (R * temp_K))) * 24 / 1000 * VS_CH4) + slurry_prod_rate * conc_fresh[['VSd_A']],
       VSnd_A = slurry_prod_rate * conc_fresh[['VSnd_A']],
       CH4_A_emis_cum = VSd_A * (exp(lnA[['VSd_A']] - E_CH4[['VSd_A']] / (R * temp_K))) * 24 / 1000,
       NH3_emis_cum = NH3_emis_rate_pit + NH3_emis_rate_floor,
       N2O_emis_cum = N2O_emis_rate,
       CH4_emis_cum = sum(rut[i_meth]) / COD_conv[['CH4']],
       CO2_emis_cum = CO2_ferm_meth_sr + respiration / COD_conv[['CO2_aer']] + alpha[['urea']] * urea / COD_conv[['CO2_ureo']],
       COD_conv_cum = sum(rut[i_meth]) + respiration + sum(rutsr),
       COD_conv_cum_meth = sum(rut[i_meth]),
       COD_conv_cum_respir = respiration,
       COD_conv_cum_sr = rutsr
     )

    return(list(derivatives, c(H2S_emis_rate = H2S_emis_rate, NH3_emis_rate_pit = NH3_emis_rate_pit,
                               NH3_emis_rate_floor = NH3_emis_rate_floor,
                               qhat = qhat, alpha = alpha, CO2_ferm = CO2_ferm, CO2_ferm_meth_sr = CO2_ferm_meth_sr, CO2_resp = respiration * COD_conv[['CO2_aer']],
                               H2S_inhib = H2S_inhib, NH3_inhib = NH3_inhib, NH4_inhib = NH4_inhib, HAC_inhib = HAC_inhib, cum_inhib = cum_inhib, 
                               rut = rut, t_run = t_run, t_batch = t_batch, conc_fresh = conc_fresh, xa_init = xa_init, 
                               xa_fresh = xa_fresh * scale['xa_fresh'], area = area, t_batch = t_batch, slurry_prod_rate = slurry_prod_rate,
                               respiration = respiration, rain = rain, evap = evap)))
     })
}

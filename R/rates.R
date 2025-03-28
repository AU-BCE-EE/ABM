rates <- function(t, y, parms, temp_C_fun = temp_C_fun, pH_fun = pH_fun, 
                  SO4_inhibition_fun = SO4_inhibition_fun, 
                  conc_fresh_fun = conc_fresh_fun, xa_fresh_fun = xa_fresh_fun) {

    y[y < 1E-10] <- 1E-10

    # need to remove slurry mass from parms to not overwrite y['slurry_mass']
    parms$slurry_mass <- NULL
    
    # all elements in parms exported to current environment
    list2env(parms, envir = environment())

    # correct slurry production rate in periods with grazing
    if(!is.null(graze_int) & any(graze_int != 0)) {
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
    
    # calculate time of a batch
    if (!is.na(wash_int)){
      batches <- c(floor((t + t_run)/(wash_int + rest_d)))
      t_batch <- (t + t_run) - batches * (wash_int + rest_d)
      if (t_batch > wash_int) t_batch <- 0
    } else {
      t_batch <- 0
    }
    
    ######REMOVE BELOW
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
    
    xa_fresh <- if (is.data.frame(xa_fresh)) {
      sapply(seq_along(grps), function(i) {
        conc <- xa_fresh_fun[[i]](t + t_run)
        return(conc)
      })
    } else{# ########## this is not needed
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
    #### replace with GREP instead or use C++
    i_meth <- grepl('^[mp]', names(qhat_opt))
    i_sr <- grepl('^sr', names(qhat_opt))
    n_mic <- length(qhat_opt)
    
    # Extract state variable values from y argument
    xa <- y[1:n_mic]
    y <- as.list(y[-c(1:n_mic)])

    # Move elements of y into rates environment
    list2env(y, envir = environment())

    # Hard-wired equilibrium constants
    log_ka <- c(NH3 = - 0.09046 - 2729.31/temp_K, 
                H2S = - 3448.7/temp_K + 47.479 - 7.5227* log(temp_K),
                HAC = -4.8288 + 21.42/temp_K)
    
    kH_oxygen <- 0.0013 * exp(1700 * ((1 / temp_K) - (1 / temp_standard))) * 32 * 1000
   
    # Hard-wire NH4+ activity coefficient
    ####### MOVE OUT OF RATES
    g_NH4 <- 0.7

    # Hydrolysis rate with Arrhenius function cpp. 
    alpha <-  Arrh_func_cpp(A, E, R, temp_K)
    names(alpha) <- names(A)
    
    if(conc_fresh$VSd <= 10E-9){
      alpha_opt_scale_type <- scale_alpha_opt[['notVSd']]
      alpha_opt_scale_CPs <- scale_alpha_opt[['CP']]
      alpha_opt_scale_CPf <- scale_alpha_opt[['CP']]
      
    } else if(conc_fresh$VSd > 10E-9){
      alpha_opt_scale_type  <- scale_alpha_opt[['VSd']]
      alpha_opt_scale_CPs <- 1
      alpha_opt_scale_CPf <- 1
    }
    
    alpha[names(alpha) != 'urea'] <- scale[['alpha_opt']] * alpha_opt_scale_type * alpha[names(alpha) != 'urea']
    alpha['CPs'] <- alpha_opt_scale_CPs * alpha['CPs']
    alpha['CPf'] <- alpha_opt_scale_CPf * alpha['CPf']
    
    # Microbial substrate utilization rate (vectorized calculation)
    qhat <- scale[['qhat_opt']] * CTM_cpp(temp_K, T_opt, T_min, T_max, qhat_opt)
    names(qhat) <- names(qhat_opt)
    
    # Decay of all microorganisms follow same kinetics with faster decay at higher temp up until 313, at which constant decay of 0.02
    decay_rate <- CTM_cpp(temp_K, 313, 273, 325, decay_rate)
    if(temp_K > 313) decay_rate <- 0.02

    # Ks temperature dependence
    ks <- ks_coefficient * (0.8157 * exp(-0.063 * temp_C)) 
    
    # NTS: Move this all out to a speciation function???
    # Chemical speciation (in rates() because is pH dependent)
    # rough approximation from Rotz et al. IFSM 2012., "markfoged" thesis, "Petersen et al. 2014", "bilds?e et al. not published", "Elzing & Aarnik 1998", and own measurements..

    HAC_frac <- 1-(1/(1 + 10^(-log_ka[['HAC']] - pH)))
    NH3_frac <- ((1/(1 + 10^(- log_ka[['NH3']] + log10(g_NH4) - pH))))
    pH_floor <- 7 # OR MOVE OUT OF RATES
    NH3_frac_floor <- ((1/(1 + 10^(- log_ka[['NH3']] + log10(g_NH4) - pH_floor))))
    
    H2S_frac <- 1 - (1/(1 + 10^(- log_ka[['H2S']] - pH))) # H2S fraction of total sulfide

    # NTS: or just add H2CO3* here?
    # NTS: need TIC production too
    
    # Using pH_LL here just to get groups
    pH_inhib  <- 0 * pH_LL + 1 
    NH3_inhib <- 0 * pH_LL + 1
    NH4_inhib <- 0 * pH_LL + 1
    HAC_inhib <- 0 * pH_LL + 1
    H2S_inhib <- 0 * pH_LL + 1


    # if pH_inhibition should be used, NH4 and NH3 inhibition is ignored and pH inhibition is used instead. 
    # inhibition is different for the microbial groups IF the inihibiton constants for are different in the grp_pars argument. 
    # Therefore the calculations are vectorized.
    ### logical expression could be faster with CPP
    if(pH_inhib_overrule){
            pH_inhib <- (1 + 2*10^(0.5* (pH_LL - pH_UL)))/(1+ 10^(pH - pH_UL) + 10^(pH_LL - pH))
    } else{
      # NH3 and NH4 inhibition      
      # Where set to 1, is still a vector to keep microbial groups in output (this is why 0 * ... bit is used)
      NH3_inhib <- ifelse(NH3_frac * TAN/slurry_mass <= ki_NH3_min, 1, exp(-2.77259 * ((NH3_frac * (TAN/(slurry_mass)) - ki_NH3_min)/(ki_NH3_max - ki_NH3_min))^2))
      NH4_inhib <- ifelse((1 - NH3_frac) * (TAN/slurry_mass) <= ki_NH4_min, 1, exp(-2.77259*(((1 - NH3_frac) * (TAN/(slurry_mass)) - ki_NH4_min)/(ki_NH4_max - ki_NH4_min))^2))
      # HAC inhibition
      # Set to 1 (no inhibition) but still a vector (below 0.05 g /kg)
      if (HAC_frac * VFA/slurry_mass >= 0.05) {
        HAC_inhib <- (2-ki_HAC/(ki_HAC+0.05)) * ki_HAC/(ki_HAC + (HAC_frac * VFA/slurry_mass))
      } else {
        HAC_inhib <- 0 * ki_HAC + 1
      }
  
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
    NH3_emis_rate_pit <- kl[['NH3']] * area * ((NH3_frac * TAN)/(slurry_mass)) * 1000 / H.NH3 # multiply by 1000 to get from g/kg to g/m3
    
    # N2O emission g(N) pr day
    N2O_emis_rate <- as.numeric(area * EF_N2O) # 
    
    # H2S emission g(S) pr day
    H2S_emis_rate <- kl[['H2S']] * area * ((H2S_frac * sulfide)/(slurry_mass)) * 1000 # multiply by 1000 to get from g/kg to g/m3

    # Respiration gCOD/d, second-order reaction where kl applies for substrate concentration of 100 g COD / kg slurry
    respiration <- 0
    sub_resp <- Cfat + CPs + CPf + RFd + starch + VSd 

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
    
    ### stoich is slow. MOVE to CPP!
    ferm <- stoich(alpha, y, conc_fresh, sub_resp, respiration)
    CO2_ferm_meth_sr <- ferm$ferm[['CO2']] * 44.01 + sum(rut[i_meth])/ferm$COD_conv_meth_CO2[[1]] + sum(rutsr)/ferm$COD_conv_sr_CO2[[1]]
    CO2_ferm <- ferm$ferm[['CO2']] * 44.01 
    
    # Derivatives, all in gCOD/d except slurry_mass = kg/d, N and S are gN or gS, VSd_A and VSnd_A are g/d
    # NTS: Some of these repeated calculations could be moved up
    # need to implement xa_bac to keep mass balance here: 
    # We need to include COD_load_rate here and possibly add COD_load_cum as state variable to enable instant output, rather than average. 
    # See around line 407 in abm()

    derivatives <- c(
       xa = scale[['yield']] * yield * rut + scale[['xa_fresh']] * xa_fresh * slurry_prod_rate - decay_rate * xa, # expands to multiple elements with element for each mic group
       slurry_mass = slurry_prod_rate + (rain - evap) * area,
       xa_aer = slurry_prod_rate * conc_fresh[['xa_aer']] + ferm[['xa_aer_rate']] - decay_rate * xa_aer,
       xa_bac = slurry_prod_rate * conc_fresh[['xa_bac']] + ferm[["xa_bac_rate"]] - decay_rate * xa_bac, # xa_bac_rate is the growth pr day (calculated in stoich function, therefore not necessary to multiply with a yield coeff)
       xa_dead = slurry_prod_rate * conc_fresh[['xa_dead']] - alpha[['xa_dead']] * xa_dead + sum(decay_rate * xa) + decay_rate * (xa_bac + xa_aer) * (1 - COD_conv[['frac_CP_xa']]), # the bacteria and methanogen mass excluding nitrogen
       RFd = slurry_prod_rate * conc_fresh[['RFd']] - alpha[['RFd']] * RFd - respiration * RFd/sub_resp,
       iNDF = slurry_prod_rate * conc_fresh[['iNDF']],
       ash = slurry_prod_rate * conc_fresh[['ash']],
       VSd = slurry_prod_rate * conc_fresh[['VSd']] - alpha[['VSd']] * VSd - respiration * VSd/sub_resp,
       starch = slurry_prod_rate * conc_fresh[['starch']] - alpha[['starch']] * starch - respiration * starch/sub_resp,
       CPs = slurry_prod_rate * conc_fresh[['CPs']] - alpha[['CPs']] * CPs - respiration * CPs/sub_resp + 
         decay_rate_xa * (xa_bac + xa_aer) * COD_conv[['frac_CP_xa']],
       CPf = slurry_prod_rate * conc_fresh[['CPf']] - alpha[['CPf']] * CPf - respiration * CPf/sub_resp,
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
       COD_conv_cum_sr = rutsr,
       COD_load_cum = slurry_prod_rate * sum(as.numeric(conc_fresh[c('starch', 'VFA', 'xa_aer', 'xa_bac', 'xa_dead', 'Cfat', 'CPs', 'CPf', 'RFd', 'iNDF', 'VSd')])) + slurry_prod_rate * sum(xa_fresh * scale[['xa_fresh']]),
       C_load_cum = slurry_prod_rate * sum(as.numeric(conc_fresh[c('starch', 'VFA', 'xa_aer', 'xa_bac', 'xa_dead', 'Cfat', 'CPs', 'CPf', 'RFd', 'iNDF', 'VSd', 'urea')]) / COD_conv[paste0('C_', c('starch', 'VFA', 'xa_aer', 'xa_bac', 'xa_dead', 'Cfat', 'CP', 'CP', 'RFd', 'iNDF', 'VSd', 'N_urea'))]) + 
         slurry_prod_rate * sum(xa_fresh * scale[['xa_fresh']] / COD_conv[['C_xa_bac']]),
       N_load_cum = slurry_prod_rate * sum(as.numeric(conc_fresh[c('CPs', 'CPf', 'TAN', 'urea')]) / c(COD_conv[['CP_N']], COD_conv[['CP_N']], 1, 1)),
       slurry_load_cum = slurry_prod_rate
     )

    return(list(derivatives, c(COD_load_rate = derivatives[['COD_load_cum']], C_load_rate = derivatives[['C_load_cum']], N_load_rate = derivatives[['N_load_cum']],
                               CH4_emis_rate = derivatives[['CH4_emis_cum']], CO2_emis_rate = derivatives[['CO2_emis_cum']], 
                               H2S_emis_rate = H2S_emis_rate, NH3_emis_rate_pit = NH3_emis_rate_pit, NH3_emis_rate_floor = NH3_emis_rate_floor,
                               N2O_emis_rate = N2O_emis_rate, CH4_A_emis_rate = derivatives[['CH4_A_emis_cum']],
                               qhat = qhat, alpha = alpha, CO2_ferm = CO2_ferm, CO2_ferm_meth_sr = CO2_ferm_meth_sr, CO2_resp = ferm[['CO2_resp']],
                               H2S_inhib = H2S_inhib, NH3_inhib = NH3_inhib, NH4_inhib = NH4_inhib,
                               TAN_min_resp = ferm[['TAN_min_resp']], TAN_min_ferm = ferm[['TAN_min_ferm']], HAC_inhib = HAC_inhib, cum_inhib = cum_inhib, 
                               rut = rut, rut_urea = rut_urea, t_run = t_run, t_batch = t_batch, conc_fresh = conc_fresh, xa_init = xa_init, 
                               xa_fresh = xa_fresh * scale[['xa_fresh']], area = area, t_batch = t_batch, slurry_prod_rate = slurry_prod_rate,
                               respiration = respiration, rain = rain, evap = evap)))
}

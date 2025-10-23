hard_pars <- function(pars){
  
  pars$g_NH4 <- 0.7
  pars$temp_standard <- 298
  pars$temp_zero <- 273
  pars$pH_floor <- 7
  pars$i_meth <- grep('^[mp]', names(pars$qhat_opt))-1 # -1 for match with CPP indexing
  pars$i_sr <- grep('^sr', names(pars$qhat_opt))-1 # -1 for match with CPP indexing
  pars$n_mic <- length(pars$qhat_opt)

  if(any(pars$conc_fresh$VSd <= 10E-9)){
    pars$alpha_opt_scale_type <- pars$scale_alpha_opt[['notVSd']]
    pars$alpha_opt_scale_CP <- pars$scale_alpha_opt[['CP']]

    
  } else if(any(pars$conc_fresh$VSd > 10E-9)){
    pars$alpha_opt_scale_type  <- pars$scale_alpha_opt[['VSd']]
    pars$alpha_opt_scale_CP <- pars$scale_alpha_opt[['CP']]
  }
  
  pars$pH_inhib  <- 0 * pars$pH_LL + 1 
  pars$NH3_inhib <- 0 * pars$pH_LL + 1
  pars$NH4_inhib <- 0 * pars$pH_LL + 1
  pars$HAC_inhib <- 0 * pars$pH_LL + 1
  pars$H2S_inhib <- 0 * pars$pH_LL + 1
  
  pars$kl[['NH3']] <- pars$kl[['NH3']] * pars$EF_NH3 

  # N2O emission g(N) pr day
  pars$N2O_emis_rate <- as.numeric(pars$area * pars$EF_N2O)
  
  # calculations moved from stoich to here to speed up model
  # with cell synthesis
  #pars$carb <- c(C6H10O5 = -1, C51H98O2 = 0, C4H6.1O1.2N = 0, 
  #               NH3 = -0.1400892, H2O = -2.3396911, C5H7O2N = 0.1400892,
  #               C2H4O2 = 1.7509058, H2 = 3.5198056, CO2 = 1.7603468)
  # without cell synthesis
  pars$carb <- c(C6H10O5 = -1, C51H98O2 = 0, C4H6.1O1.2N = 0, 
                 NH3 = 0, H2O = -2.992714, C5H7O2N = 0,
                 C2H4O2 = 2.003743, H2 = 3.985329, CO2 = 1.993114)
  
  
  # with cell synthesis
  #pars$pro <- c(C6H10O5 = 0, C51H98O2 = 0, C4H6.1O1.2N = -1, 
  #              NH3 = 0.9098195, H2O = -2.91011939, C5H7O2N = 0.09091805,
  #              C2H4O2 = 1.32238963, H2 = 1.24392800, CO2 = 0.53512875)
  # without cell synthesis
  pars$pro <- c(C6H10O5 = 0, C51H98O2 = 0, C4H6.1O1.2N = -1, 
                NH3 = 1, H2O = -3.2521987, C5H7O2N = 0,
                C2H4O2 = 1.4774587, H2 = 1.402611, CO2 = 0.6004123)
  # with cell synthesis
  #pars$lip <- c(C6H10O5 = 0, C51H98O2 = -1, C4H6.1O1.2N = 0, 
  #              NH3 = -0.9111515, H2O = -41.2980881, C5H7O2N = 0.9111515,
  #              C2H4O2 = 23.7180355, H2 = 40.9726177, CO2 = 1.4659059)
  # without cell synthesis
  pars$lip <- c(C6H10O5 = 0, C51H98O2 = -1, C4H6.1O1.2N = 0, 
                NH3 = 0, H2O = -45.4119381, C5H7O2N = 0,
                C2H4O2 = 25.3171914, H2 = 43.7334158, CO2 = 0.3667798)
  
  # VSd without cell synthesis (based on C15.89H26.3O9.2N)  
  pars$OM <- c(C6H10O5 = 0, C51H98O2 = 0, C4H6.1O1.2N = 0, 
                NH3 = 1, H2O = -11.284, C5H7O2N = 0,
                C2H4O2 = 5.64795, H2 = 11.638, CO2 = 4.5941)
  
  # with cell synthesis
  #carb_resp <- c(C6H12O6 = -1, C51H98O2 = 0, C4H6.1O1.2N = 0, 
  #                    NH4 = -0.78, HCO3 = -0.78, O2 =  -2.1, H2O = 5.22, C5H7O2N = 0.78,
  #                    CO2 = 2.88) # 
  carb_resp <- c(C6H12O6 = -1, C51H98O2 = 0, C4H6.1O1.2N = 0, 
                      NH4 = 0, HCO3 = 0, O2 = -6, H2O = 6, C5H7O2N = 0,
                      CO2 = 6) # 
  
  pro_resp <- c(C6H12O6 = 0, C51H98O2 = 0, C4H6.1O1.2N = -1, 
                     NH4 = 1, HCO3 = 0, O2 = -4.175, H2O = 1.55, C5H7O2N = 0,
                     CO2 = 4) #
  
  lip_resp <- c(C6H12O6 = 0, C51H98O2 = -1, C4H6.1O1.2N = 0, 
                     NH4 = 0, HCO3 = 0, O2 = -74.5, H2O = 49, C5H7O2N = 0,
                     CO2 = 51) #
  ## THIS IS WHERE I AM
  OM_resp <- c(C6H12O6 = 0, C51H98O2 = -1, C4H6.1O1.2N = 0, 
                NH4 = 1, HCO3 = 0, O2 = -17.865, H2O = 13.15, C5H7O2N = 0,
                CO2 = 15.89) #

  pars$carb_resp <- carb_resp # 
  pars$pro_resp <- pro_resp # 
  pars$lip_resp <- lip_resp #
  pars$OM_resp <- OM_resp # 
  
  pars$ace <- c(H2 = 0, C2H4O2 = -1, CO2 = 1, CH4 = 1, H2O = 0) # moles stoichiometry 
  pars$hyd <- c(H2 = -1, C2H4O2 = 0, CO2 = -1/4, CH4 = 1/4, H2O = 2/4)
  pars$ace_sr <- c(H2 = 0, C2H4O2 = -1, H2SO4 = -1, CO2 = 2, H2O = 2, H2S = 1) 
  pars$hyd_sr <- c(H2 = -1, C2H4O2 = 0, H2SO4 = -1/4, CO2 = 0, H2O = 1, H2S = 1/4) 
  


  return(pars)
}
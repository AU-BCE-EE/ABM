hard_pars <- function(pars){
  
  pars$g_NH4 <- 0.7
  pars$temp_standard <- 298
  pars$temp_zero <- 273
  pars$pH_floor <- 7
  pars$i_meth <- grep('^[mp]', names(pars$qhat_opt))-1 # -1 for match with CPP indexing
  pars$i_sr <- grep('^sr', names(pars$qhat_opt))-1 # -1 for match with CPP indexing
  pars$n_mic <- length(pars$qhat_opt)
  pars$rut <- NA * pars$qhat_opt
  
  if(pars$conc_fresh$VSd <= 10E-9){
    pars$alpha_opt_scale_type <- pars$scale_alpha_opt[['notVSd']]
    pars$alpha_opt_scale_CP <- pars$scale_alpha_opt[['CP']]

    
  } else if(pars$conc_fresh$VSd > 10E-9){
    pars$alpha_opt_scale_type  <- pars$scale_alpha_opt[['VSd']]
    pars$alpha_opt_scale_CP <- 1
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
  pars$carb <- c(C6H10O5 = -1, C51H98O2 = 0, C4H6.1O1.2N = 0, 
                 NH3 = -0.1400892, H2O = -2.3396911, C5H7O2N = 0.1400892,
                 C2H4O2 = 1.7509058, H2 = 3.5198056, CO2 = 1.7603468)
  
  pars$pro <- c(C6H10O5 = 0, C51H98O2 = 0, C4H6.1O1.2N = -1, 
                NH3 = 0.9098195, H2O = -2.91011939, C5H7O2N = 0.09091805,
                C2H4O2 = 1.32238963, H2 = 1.24392800, CO2 = 0.53512875)
  
  pars$lip <- c(C6H10O5 = 0, C51H98O2 = -1, C4H6.1O1.2N = 0, 
                NH3 = -0.9111515, H2O = -41.2980881, C5H7O2N = 0.9111515,
                C2H4O2 = 23.7180355, H2 = 40.9726177, CO2 = 1.4659059)
  
  pars$carb_resp <- c(C6H12O6 = -1, C51H98O2 = 0, C4H6.1O1.2N = 0, 
                      NH4 = -0.78, HCO3 = -0.78, O2 =  -2.1, H2O = 5.22, C5H7O2N = 0.78,
                      CO2 = 2.88)
  pars$pro_resp <- c(C6H12O6 = 0, C51H98O2 = 0, C4H6.1O1.2N = -1, 
                     NH4 = 0.45725, HCO3 = 0.45725, O2 = -1.46125, H2O = 0.007249637, C5H7O2N = 0.54275,
                     CO2 = 0.829)
  pars$lip_resp <- c(C6H12O6 = 0, C51H98O2 = -1, C4H6.1O1.2N = 0, 
                     NH4 = -9.425, HCO3 = -9.425, O2 = -25.375, H2O = 39.575, C5H7O2N = 9.425,
                     CO2 = 13.3) 
  
  pars$ace <- c(H2 = 0, C2H4O2 = -1, CO2 = 1, CH4 = 1, H2O = 0)
  pars$hyd <- c(H2 = -1, C2H4O2 = 0, CO2 = -1/4, CH4 = 1/4, H2O = 2/4)
  pars$ace_sr <- c(H2 = 0, C2H4O2 = -1, H2SO4 = -1, CO2 = 2, H2O = 2, H2S = 1) 
  pars$hyd_sr <- c(H2 = -1, C2H4O2 = 0, H2SO4 = -1/4, CO2 = 0, H2O = 1, H2S = 1/4) 
  
  return(pars)
}
resp_stoich <- function(conc_fresh, carb_resp, pro_resp, lip_resp){
  
  # mole fermented per day: coefficients in the end has unit of mole/gCOD and "mol.carb, mol.pro, mol.lip" is mole /day
  
  sub_resp <- conc_fresh$starch + conc_fresh$Cfat + conc_fresh$RFd + conc_fresh$CPf + conc_fresh$CPs + conc_fresh$VSd
  
  if (conc_fresh['VSd'] <= 1e-10) { # macro molecular approach
    mol.carb_resp <- (conc_fresh$RFd + conc_fresh$starch) / sub_resp * 0.005208333 * carb_resp
    mol.pro_resp  <- (conc_fresh$CPs + conc_fresh$CPf) / sub_resp * 0.00748503 * pro_resp
    mol.lip_resp  <- conc_fresh$Cfat / sub_resp * 0.0004194631 * lip_resp
  } else { # VSd approach
    mol.carb_resp <- (conc_fresh$VSd) / sub_resp * 0.002753327 * carb_resp
    mol.pro_resp  <- (conc_fresh$VSd)/ sub_resp * 0.00176104 * pro_resp
    mol.lip_resp  <- (conc_fresh$VSd)/ sub_resp * 9.902938e-05 * lip_resp
  }
  
  resp <- mol.carb_resp + mol.pro_resp + mol.lip_resp
  xa_aer_rate <- resp[['C5H7O2N']] * 160 # need to be multiplied by respiration in rates (gCOD/day), 113.113 * 1.414515  # aerobic bacteria growth in gCOD/day, 113.113 is g pr mol biomass, and 1.41.. is gCOD / g biomass 
  TAN_min_resp <- resp[['NH4']] * 14.007 # need to be multiplied by respiration in rates (gCOD/day), g N / day
  CO2_resp <- (resp[['CO2']] + resp[['HCO3']]) * 44.01 # gCO2/gCOD, need to be multiplied by respiration in rates (gCOD/day)

  return(list(CO2_resp = CO2_resp, xa_aer_rate = xa_aer_rate, TAN_min_resp = TAN_min_resp))
}
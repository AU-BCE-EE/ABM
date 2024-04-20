stoich <- function(alpha, y, conc_fresh, sub_resp, respiration){

# mole fermented per day: coefficients in the end has unit of mole/gCOD and "mol.carb, mol.pro, mol.lip" is mole/day
mol.carb <- (alpha['RFd'] * y$RFd + alpha['starch'] * y$starch) * 0.005208333
mol.pro <- (alpha['CP'] * y$CP) * 0.00748503
mol.lip <- (alpha['Cfat'] * y$Cfat) * 0.0004194631

mol.carb_resp <- respiration * (y$RFd + y$starch)/sub_resp * 0.005208333
mol.pro_resp <- respiration * y$CP/sub_resp * 0.00748503
mol.lip_resp <- respiration * y$Cfat/sub_resp * 0.0004194631

# if VSd is used. based on composition of degradable cattle manure excluding vfa (Appendix 1, ABM paper)
if(conc_fresh[['VSd']] > 1e-10){
  mol.carb <- alpha['VSd'] * y$VSd * 0.002753327
  mol.pro <- alpha['VSd'] * y$VSd * 0.00176104
  mol.lip <- alpha['VSd'] * y$VSd * 9.902938e-05
  mol.carb_resp <- respiration * y$VSd * 0.002753327
  mol.pro_resp <- respiration * y$VSd * 0.00176104
  mol.carb_resp <- respiration * y$VSd * 9.902938e-05
}

# stoichiometry. Now assuming cell synthesis, mole/day
carb <- c(C6H10O5 = -1, C51H98O2 = 0, C4H6.1O1.2N = 0, 
         NH3 = -0.1400892, H2O = -2.3396911, C5H7O2N = 0.1400892,
         C2H4O2 = 1.7509058, H2 = 3.5198056, CO2 = 1.7603468) * mol.carb

# stoichiometry. Now assuming cell synthesis, mole/day
pro <- c(C6H10O5 = 0, C51H98O2 = 0, C4H6.1O1.2N = -1, 
         NH3 = 0.9098195, H2O = -2.91011939, C5H7O2N = 0.09091805,
         C2H4O2 = 1.32238963, H2 = 1.24392800, CO2 = 0.53512875) * mol.pro

# stoichiometry. Now assuming cell synthesis, mole/day
lip <- c(C6H10O5 = 0, C51H98O2 = -1, C4H6.1O1.2N = 0, 
         NH3 = -0.9111515, H2O = -41.2980881, C5H7O2N = 0.9111515,
         C2H4O2 = 23.7180355, H2 = 40.9726177, CO2 = 1.4659059) * mol.lip

ferm <- carb + pro + lip

ace <- c(H2 = 0, C2H4O2 = -1, CO2 = 1, CH4 = 1, H2O = 0) * ferm['C2H4O2']  
hyd <- c(H2 = -1, C2H4O2 = 0, CO2 = -1/4, CH4 = 1/4, H2O = 2/4) * ferm['H2']      

ace_sr <- c(H2 = 0, C2H4O2 = -1, H2SO4 = -1, CO2 = 2, H2O = 2, H2S = 1) * ferm['C2H4O2'] 
hyd_sr <- c(H2 = -1, C2H4O2 = 0, H2SO4 = -1/4, CO2 = 0, H2O = 1, H2S = 1/4) * ferm['H2'] 

# combine methanogenesis stoichiometry for calculating CO2 conv factor below. 
meth <- c(C2H4O2 = ace[['C2H4O2']], H2 = hyd[['H2']], 
         CH4 = ace[['CH4']] + hyd[['CH4']],
         CO2 = hyd[['CO2']] + ace[['CO2']])

# combine sulfate reduction stoichiometry for calculating CO2 conv factor below. 
sr <-  c(C2H4O2 = ace_sr[['C2H4O2']], H2 = hyd_sr[['H2']], 
         H2SO4 = ace_sr[['H2SO4']] + hyd_sr[['H2SO4']],
         H2S = ace_sr[['H2S']] + hyd_sr[['H2S']],
         CO2 = ace_sr[['CO2']])

xa_bac_rate <- ferm[['C5H7O2N']] * 113.113 * 1.414515  # fermentative and hydrolytic bacteria growth in gCOD/day, 113.113 is g pr mol biomass, and 1.41.. is gCOD / g biomass 
VFA_H2 <- ferm[['C2H4O2']] * 60.052 * 1.065743 + ferm[['H2']] * 2.016 * 7.936508 # acetate and H2 COD from hydrolysis + fermentation reactions.
TAN_min_ferm <- ferm[['NH3']] * 14.007 # g N-NH3

# COD_conversion from moles to gCOD (acetate and H2) or moles to g (for CO2), 
# coefficients in gCOD/mol on reaction side and g/mol for CO2. 
# conversion factor (COD_conv_meth_CO2) has unit gCOD consumed/gCO2 produced
# conversion factor (COD_conv_sr_CO2) has unit gCOD consumed/gCO2 produced

COD_conv_meth_CO2 <- -(meth['C2H4O2'] * 64 + meth['H2'] * 16)/(meth['CO2'] * 44.01)
COD_conv_sr_CO2 <- -(sr['C2H4O2'] * 64 + sr['H2'] * 16)/(sr['CO2'] * 44.01)

# stoichiometry  for respiration. Now assuming cell synthesis, mole/day
# energy fraction going to growth, fs0, is 0.65 for aerobic bacteria, Rittman, it is build into the equations below 
carb_resp <- c(C6H12O6 = -1, C51H98O2 = 0, C4H6.1O1.2N = 0, 
          NH4 = -0.78, HCO3 = -0.78, O2 =  -2.1, H2O = 5.22, C5H7O2N = 0.78,
          CO2 = 2.88) * mol.carb_resp

# stoichiometry for respiration. Now assuming cell synthesis, mole/day
pro_resp <- c(C6H12O6 = 0, C51H98O2 = 0, C4H6.1O1.2N = -1, 
         NH4 = 0.45725, HCO3 = 0.45725, O2 = -1.46125, H2O = 0.007249637, C5H7O2N = 0.54275,
         CO2 = 0.829) * mol.pro_resp

# stoichiometry  for respiration. Now assuming cell synthesis, mole/day
lip_resp <- c(C6H12O6 = 0, C51H98O2 = -1, C4H6.1O1.2N = 0, 
         NH4 = -9.425, HCO3 = -9.425, O2 = -25.375, H2O = 39.575, C5H7O2N = 9.425,
         CO2 = 13.3) * mol.lip_resp

resp <- carb_resp + pro_resp + lip_resp
xa_aer_rate <- resp[['C5H7O2N']] * 113.113 * 1.414515  # aerobic bacteria growth in gCOD/day, 113.113 is g pr mol biomass, and 1.41.. is gCOD / g biomass 
TAN_min_resp <- resp[['NH4']] * 14.007 # g N / day
CO2_resp <- (resp[['CO2']] + resp[['HCO3']]) * 44.01 # g CO2 / day

return(list(ferm = ferm, COD_conv_meth_CO2 = COD_conv_meth_CO2, COD_conv_sr_CO2 = COD_conv_sr_CO2, 
            xa_bac_rate = xa_bac_rate, xa_aer_rate = xa_aer_rate, 
            TAN_min_ferm = TAN_min_ferm, TAN_min_resp = TAN_min_resp,
            VFA_H2 = VFA_H2, CO2_resp = CO2_resp))

}


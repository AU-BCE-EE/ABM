stoich <- function(alpha, y, conc_fresh){

# mole fermented per day: coefficient has unit of mole/gCOD
mol.carb <- (alpha['RFd'] * y$RFd + alpha['starch'] * y$starch) * 0.005208333
mol.pro <- (alpha['CP'] * y$CP) * 0.00748503
mol.lip <- (alpha['Cfat'] * y$Cfat) * 0.0004194631

# if VSd is used. based on composition of degradable cattle manure excluding vfa (Appendix 1, ABM paper)
if(conc_fresh[['VSd']] > 1e-10){
  mol.carb <- alpha['VSd'] * y$VSd * 0.002753327
  mol.pro <- alpha['VSd'] * y$VSd * 0.00176104
  mol.lip <- alpha['VSd'] * y$VSd * 9.902938e-05
}

# stoichiometry. Assuming no cell synthesis
carb <- c(C6H10O5 = -1, C51H98O2 = 0, C4H6.1O1.2N = 0, 
         NH3 = 0, H2O = -3, C5H7O2N = 0,
         C2H4O2 = 2, H2 = 4, CO2 = 2) * mol.carb

# stoichiometry. Assuming no cell synthesis
pro <- c(C6H10O5 = 0, C51H98O2 = 0, C4H6.1O1.2N = -1, 
         NH3 = 1, H2O = -3.549196, C5H7O2N = 0,
         C2H4O2 = 1.6252168, H2 = 1.8484773, CO2 = 0.7492813) * mol.pro

# stoichiometry. Assuming no cell synthesis
lip <- c(C6H10O5 = 0, C51H98O2 = -1, C4H6.1O1.2N = 0, 
         NH3 = 0, H2O = -45.4119381, C5H7O2N = 0,
         C2H4O2 = 25.3171914, H2 = 43.7334158, CO2 = 0.3667798) * mol.lip

ferm <- carb + pro + lip

ace <- c(H2 = 0, C2H4O2 = -1, CO2 = 1, CH4 = 1, H2O = 0) * ferm['C2H4O2']  
hyd <- c(H2 = -1, C2H4O2 = 0, CO2 = -1/4, CH4 = 1/4, H2O = 2/4) * ferm['H2']      

ace_sr <- c(H2 = 0, C2H4O2 = -1, H2SO4 = -1, CO2 = 2, H2O = 2, H2S = 1) * ferm['C2H4O2'] 
hyd_sr <- c(H2 = -1, C2H4O2 = 0, H2SO4 = -1/4, CO2 = 0, H2O = 4/4, H2S = 1) * ferm['H2'] 

# combine methanogenesis stoichiometry for calculating CO2 conv factor below. 
meth <- c(C2H4O2 = ace[['C2H4O2']], H2 = hyd[['H2']], 
         CH4 = ace[['CH4']] + hyd[['CH4']],
         CO2 = hyd[['CO2']] + ace[['CO2']])

# combine sulfate reduction stoichiometry for calculating CO2 conv factor below. 
sr <-  c(C2H4O2 = ace_sr[['C2H4O2']], H2 = hyd_sr[['H2']], 
         H2SO4 = ace_sr[['H2SO4']] + hyd_sr[['H2SO4']],
         H2S = ace_sr[['H2S']] + hyd_sr[['H2S']],
         CO2 = ace_sr[['CO2']])

# COD_conversion from moles to gCOD (acetate and H2) or moles to g (for CO2), 
# coefficients in gCOD/mol on reaction side and g/mol for CO2. 
# conversion factor (COD_conv_meth_CO2) has unit gCOD consumed/gCO2 produced
# conversion factor (COD_conv_sr_CO2) has unit gCOD consumed/gCO2 produced

COD_conv_meth_CO2 <- -(meth['C2H4O2'] * 64 + meth['H2'] * 16)/(meth['CO2'] * 44.01)
COD_conv_sr_CO2 <- -(sr['C2H4O2'] * 64 + sr['H2'] * 16)/(sr['CO2'] * 44.01)

return(list(ferm = ferm, COD_conv_meth_CO2 = COD_conv_meth_CO2, COD_conv_sr_CO2 = COD_conv_sr_CO2))

}


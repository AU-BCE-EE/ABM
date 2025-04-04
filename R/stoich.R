stoich <- function(alpha, y, conc_fresh, sub_resp, respiration, 
                   carb, pro, lip, carb_resp, pro_resp, lip_resp, ace, hyd,
                   ace_sr, hyd_sr){

# mole fermented per day: coefficients in the end has unit of mole/gCOD and "mol.carb, mol.pro, mol.lip" is mole /day
mol.carb <- (alpha['RFd'] * y['RFd'] + alpha['starch'] * y['starch']) * 0.005208333
mol.pro <- (alpha['CPs'] * y['CPs'] + alpha['CPf'] * y['CPf']) * 0.00748503
mol.lip <- (alpha['Cfat'] * y['Cfat']) * 0.0004194631

mol.carb_resp <- respiration * (y['RFd'] + y['starch'])/sub_resp * 0.005208333
mol.pro_resp <- (respiration * y['CPs']/sub_resp + respiration * y['CPf']/sub_resp)  * 0.00748503
mol.lip_resp <- respiration * y['Cfat']/sub_resp * 0.0004194631

# if VSd is used. based on composition of degradable cattle manure excluding 
# vfa (Appendix 1, ABM paper)
if(conc_fresh[['VSd']] > 1e-10){
  mol.carb <- alpha['VSd'] * y['VSd'] * 0.002753327
  mol.pro <- alpha['VSd'] * y['VSd'] * 0.00176104
  mol.lip <- alpha['VSd'] * y['VSd'] * 9.902938e-05
  mol.carb_resp <- respiration * y['VSd']/sub_resp  * 0.002753327
  mol.pro_resp <- respiration * y['VSd']/sub_resp  * 0.00176104
  mol.lip_resp <- respiration * y['VSd']/sub_resp  * 9.902938e-05
}

# stoichiometry. Now assuming cell synthesis, mole/day
carb <- carb * mol.carb

# stoichiometry. Now assuming cell synthesis, mole/day
pro <- pro * mol.pro

# stoichiometry. Now assuming cell synthesis, mole/day
lip <- lip * mol.lip

ferm <- carb + pro + lip

ace <- ace * ferm['C2H4O2']  
hyd <- hyd * ferm['H2']      

ace_sr <- ace_sr * ferm['C2H4O2'] 
hyd_sr <- hyd_sr * ferm['H2'] 

# combine methanogenesis stoichiometry for calculating CO2 conv factor below. 
meth <- c(C2H4O2 = ace[['C2H4O2']], H2 = hyd[['H2']], 
         CH4 = ace[['CH4']] + hyd[['CH4']],
         CO2 = hyd[['CO2']] + ace[['CO2']])

# combine sulfate reduction stoichiometry for calculating CO2 conv factor below. 
sr <-  c(C2H4O2 = ace_sr[['C2H4O2']], H2 = hyd_sr[['H2']], 
         H2SO4 = ace_sr[['H2SO4']] + hyd_sr[['H2SO4']],
         H2S = ace_sr[['H2S']] + hyd_sr[['H2S']],
         CO2 = ace_sr[['CO2']])

xa_bac_rate <- ferm[['C5H7O2N']] * 160 # 113.113 g/mol * 1.414515 gCOD/g # fermentative and hydrolytic bacteria growth in gCOD/day, 113.113 is g pr mol biomass, and 1.41.. is gCOD / g biomass 
VFA_H2 <- ferm[['C2H4O2']] * 64 + ferm[['H2']] * 16 # 60.052 * 1.065743 = 64, 2.016 * 7.936508 = 16, acetate and H2 COD from hydrolysis + fermentation reactions.This is energy that can be consumed by methanogenesis and sulfate reduction
TAN_min_ferm <- ferm[['NH3']] * 14.007 # g N-NH3, how much TAN is produced from fermentation (mineralization of protein)

# COD_conv: from moles to gCOD (acetate and H2) or moles to g (for CO2), 
# coefficients in gCOD/mol on reaction side and g/mol for CO2. 
# conversion factor (COD_conv_meth_CO2) has unit gCOD consumed/gCO2 produced
# conversion factor (COD_conv_sr_CO2) has unit gCOD consumed/gCO2 produced

COD_conv_meth_CO2 <- -(meth['C2H4O2'] * 64 + meth['H2'] * 16)/(meth['CO2'] * 44.01) # 64 is gCOD/mol acetate, 16 is gCOD/mol H2, 44 is gCO2/mol
COD_conv_sr_CO2 <- -(sr['C2H4O2'] * 64 + sr['H2'] * 16)/(sr['CO2'] * 44.01)

# stoichiometry  for respiration. Now assuming cell synthesis, mole/day
# energy fraction going to growth, fs0, is 0.65 for aerobic bacteria, Rittman, it is build into the equations below 
carb_resp <- carb_resp * mol.carb_resp

# stoichiometry for respiration. Now assuming cell synthesis, mole/day
pro_resp <- pro_resp * mol.pro_resp

# stoichiometry  for respiration. Now assuming cell synthesis, mole/day
lip_resp <- lip_resp * mol.lip_resp

resp <- carb_resp + pro_resp + lip_resp
xa_aer_rate <- resp[['C5H7O2N']] * 160 # 113.113 * 1.414515  # aerobic bacteria growth in gCOD/day, 113.113 is g pr mol biomass, and 1.41.. is gCOD / g biomass 
TAN_min_resp <- resp[['NH4']] * 14.007 # g N / day
CO2_resp <- (resp[['CO2']] + resp[['HCO3']]) * 44.01 # g CO2 / day

return(list(ferm = ferm, COD_conv_meth_CO2 = COD_conv_meth_CO2, COD_conv_sr_CO2 = COD_conv_sr_CO2, 
            xa_bac_rate = xa_bac_rate, xa_aer_rate = xa_aer_rate, 
            TAN_min_ferm = TAN_min_ferm, TAN_min_resp = TAN_min_resp,
            VFA_H2 = VFA_H2, CO2_resp = CO2_resp))

}


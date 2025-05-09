# Concept for stoich()
# S. Hafner & F. Dalby

# Substrate state variables
source('../../R/predFerm.R')
source('../../R/readFormula.R')

# unit of y is moles here. 
y <- c(xa = 0,
       slurry_mass = 0, 
       xa_aer = 0,
       xa_bac = 0,
       xa_dead = 0, 
       RFd = 1,
       iNDF = 0,
       ash = 0,
       VSd = 0,
       starch = 0,
       CPs = 0,
       CPf = 0,
       Cfat = 0,
       VFA = 0, 
       urea = 0, 
       TAN = 0, 
       sulfate = 0, 
       sulfide = 0,
       NH3_emis_cum = 0, 
       N2O_emis_cum = 0, 
       CH4_emis_cum = 0, 
       CO2_emis_cum = 0, 
       COD_conv_cum = 0, 
       COD_conv_cum_meth = 0, 
       COD_conv_cum_respir = 0, 
       COD_conv_cum_sr = 0,
       COD_load_cum = 0,
       C_load_cum = 0,
       N_load_cum = 0,
       slurry_load_cum = 0)

fs = 0.1175#
acefrac = 0.66
# state variables that are involved in hydrolysis/fermentation part (both consumption and production)
#y_names_alpha <- c('xa_dead','xa_bac','RFd','VSd','starch','CPs','CPf','Cfat','VFA','urea','TAN','CO2_emis_cum')
y_names <- names(y)
  
# calc molar stoichiometry for compounds being hydrolyzed through alpha

## NTS need to convert all below to gCOD or g or gN or gS directly in the stoichiometry. Then we can avoind the sweep with mole_to_gCOD further down.
## this conversion is not done yet.

# xa_dead is tricky because we have tried to separate it from Nitrogen part. All the COD of xa_dead there goes into VFA here.
xa_dead <- c(xa_dead = -1, CH3COOH = 1)

starch <- predFerm('C6H10O5', acefrac = acefrac, fs = fs)/predFerm('C6H10O5', acefrac = acefrac, fs = fs)['C6H10O5']*-1
names(starch)[names(starch) == "C6H10O5"] <- "starch"

Cfat <- predFerm('C51H98O6', acefrac = acefrac, fs = fs)/predFerm('C51H98O6', acefrac = acefrac, fs = fs)['C51H98O6']*-1
names(Cfat)[names(Cfat) == "C51H98O6"] <- "Cfat"

CPs <- predFerm('C4H6.1O1.2N', acefrac = acefrac, fs = fs)/predFerm('C4H6.1O1.2N', acefrac = acefrac, fs = fs)['C4H6.1O1.2N']*-1
names(CPs)[names(CPs) == "C4H6.1O1.2N"] <- "CPs"

CPf <- predFerm('C4H6.1O1.2N', acefrac = acefrac, fs = fs)/predFerm('C4H6.1O1.2N', acefrac = acefrac, fs = fs)['C4H6.1O1.2N']*-1
names(CPf)[names(CPf) == "C4H6.1O1.2N"] <- "CPf"

RFd <- predFerm('C6H10O5', acefrac = acefrac, fs = fs)/predFerm('C6H10O5', acefrac = acefrac, fs = fs)['C6H10O5']*-1
names(RFd)[names(RFd) == "C6H10O5"] <- "RFd"

VSd <- predFerm('C15.9H26.34O9.2N', acefrac = acefrac, fs = fs)/predFerm('C15.9H26.34O9.2N', acefrac = acefrac, fs = fs)['C15.9H26.34O9.2N']*-1
names(VSd)[names(VSd) == "C15.9H26.34O9.2N"] <- "VSd"

# urea is also different since it is just the rut urea_N so it is assumed to be converted 1 to 1 to TAN
urea <- c(urea = -1, NH3 = 1, CO2 = 0.5, H2O = -0.5) # this is per urea N

# add names to stoichiometry
# names_stoich <- c(names(RFd), 'starch', 'CPs', 'CPf','Cfat','VSd','urea', 'xa_dead')
names_stoich <- c(names(RFd), 'starch', 'CPs', 'CPf','Cfat', 'VFA', 'VSd','urea', 'xa_dead','sulfide','sulfate','xa','slurry_mass', 'xa_aer','iNDF','ash',
                  'NH3_emis_cum','N2O_emis_cum','CH4_emis_cum','COD_conv_cum','COD_conv_cum_meth','COD_conv_cum_respir','COD_conv_cum_sr', 'COD_load_cum',
                  'C_load_cum', 'N_load_cum', 'slurry_load_cum')
# make matrix with state vars in rows and stoich in cols
cc <- matrix(nrow = length(y_names), ncol = length(names_stoich), dimnames = list(y_names, names_stoich))

# fill in with stoichiometry for compounds being consumed only.
cc['RFd', names(RFd)] <- RFd
cc['starch', names(starch)] <- starch
cc['CPs', names(CPs)] <- CPs
cc['CPf', names(CPf)] <- CPf
cc['VSd', names(VSd)] <- VSd
cc['Cfat', names(Cfat)] <- Cfat
cc['xa_dead', names(xa_dead)] <- xa_dead
cc['urea', names(urea)] <- urea

# fill out all NA with 0
cc[is.na(cc)] <- 0
colnames(cc)[colnames(cc) == "NH3"] <- "TAN"
colnames(cc)[colnames(cc) == "C5H7O2N"] <- "xa_bac"
colnames(cc)[colnames(cc) == "CO2"] <- "CO2_emis_cum"

# convert moles to gCOD for the compounds involved in hydrolysis/fermentation to VFA
# C5H7O2N = 113.113 g/mol * 1.414515 gCOD/g
moles_to_gCOD <- c(H2O = 16, 
                   starch = biogas::molMass('C6H10O5') * biogas::calcCOD('C6H10O5'), 
                   RFd = biogas::molMass('C6H10O5') * biogas::calcCOD('C6H10O5'), 
                   TAN = 14.007, 
                   H2 = biogas::molMass('H2') * biogas::calcCOD('H2'), 
                   CO2_emis_cum = 44.01, 
                   CH3COOH = biogas::molMass('CH3COOH') * biogas::calcCOD('CH3COOH'), 
                   xa_bac = biogas::molMass('C5H7O2N') * biogas::calcCOD('C5H7O2N'), 
                   xa_aer = biogas::molMass('C5H7O2N') * biogas::calcCOD('C5H7O2N'),
                   xa_dead = biogas::molMass('C5H7O2N') * biogas::calcCOD('C5H7O2N'),
                   xa = biogas::molMass('C5H7O2N') * biogas::calcCOD('C5H7O2N'),
                   CPs = biogas::molMass('C4H6.1O1.2N') * biogas::calcCOD('C4H6.1O1.2N'), 
                   CPf = biogas::molMass('C4H6.1O1.2N') * biogas::calcCOD('C4H6.1O1.2N'), 
                   Cfat = biogas::molMass('C51H98O6') * biogas::calcCOD('C51H98O6'), 
                   VSd = biogas::molMass('C15.9H26.34O9.2N') * biogas::calcCOD('C15.9H26.34O9.2N'),
                   urea = 28.014)
# add names to moles_to_gCOD so it matches with colnames in cc
# then just set those elements to 0
new_names <- setdiff(colnames(cc), names(moles_to_gCOD))
moles_to_gCOD[new_names] <- 0

# reorder to ensure same order
moles_to_gCOD <- moles_to_gCOD[colnames(cc)]
# check consistency in lengths
length(moles_to_gCOD) == ncol(cc)
# change unit of matrix to gCOD/mol or gN/mol or g/mol. 
# Need to convert now because we need to merge H2 and CH3COOH and they 
# should be in COD units before doign that. 
# converting coefficients to pr gCOD
cc_gCOD_mole <- sweep(cc, 2, moles_to_gCOD, `*`) # 
cc_COD <- cc_gCOD_mole
cc_COD['xa_dead',] <- -cc_gCOD_mole['xa_dead',]/cc_gCOD_mole['xa_dead','xa_dead']
cc_COD['RFd',] <- -cc_gCOD_mole['RFd',]/cc_gCOD_mole['RFd','RFd']
cc_COD['VSd',] <- -cc_gCOD_mole['VSd',]/cc_gCOD_mole['VSd','VSd']
cc_COD['starch',] <- -cc_gCOD_mole['starch',]/cc_gCOD_mole['starch','starch']
cc_COD['CPs',] <- -cc_gCOD_mole['CPs',]/cc_gCOD_mole['CPs','CPs']
cc_COD['CPf',] <- -cc_gCOD_mole['CPf',]/cc_gCOD_mole['CPf','CPf']
cc_COD['Cfat',] <- -cc_gCOD_mole['Cfat',]/cc_gCOD_mole['Cfat','Cfat']
cc_COD['urea',] <- -cc_gCOD_mole['urea',]/cc_gCOD_mole['urea','urea']

# remove water and rename relevant column names
cc_COD <- cc_COD[, colnames(cc_COD) != 'H2O']
# H2 + CH3COOH = VFA column
cc_COD[, 'VFA'] <- cc_COD[,'H2'] + cc_COD[, 'CH3COOH']
# remove H2 and CH3COOH
cc_COD <- cc_COD[, !colnames(cc_COD) %in% c('H2','CH3COOH')]
# reorder column names in same order as rows
cc_COD <- cc_COD[, y_names]
# transpose it, such that rows = state variable and columns are coefficients multiplied with column state var in dot-product
cc_t <- t(cc_COD)
cc_t
# cc_t should be passed into rates environment
# alpha should be calculated and depends on e.g. 
# temperature, which is dynamic and has to be in rates
# creating alpha dummy. Als for compounds not being hydrolyzed. It does
# not matter because they have no coefficients anyway (see VFA or bac column)
alpha <- rep(1, length(y_names_alpha)) 
# multiplying with alpha along each row
cc_alpha <- sweep(cc_t, 1, alpha, `*`)
# this is the change in the state variables when only accounting for hydrolysis
# and since y[y_names_alpha] has units of moles here.
# Ideally should hard code the stoich in units of COD instead 
# When 1 gCOD is consumed how much is produced of the stoichiometry 
alpha_derivative <-  cc_alpha %*% y
alpha_derivative
#### Making slurry addition matrix
# now we create the fresh slurry addition matrix and set it to 0 initially
cc_fresh <- cc_t
cc_fresh[,] <- 0
# get conc_fresh this can be set outside rates_cpp if conc_fresh is constant
# otherwise it needs to be done inside rates_cpp. fresh concentration go into coefficient matrix
# and this is dotted with slurry_prod_rate
# dummy conc_fresh
conc_fresh <- ABM::man_pars2.0$conc_fresh
matching_names <- intersect(names(conc_fresh), rownames(cc_fresh))
cc_fresh[cbind(matching_names, matching_names)] <- unlist(conc_fresh[matching_names])
cc_fresh['slurry_mass', 'slurry_mass'] <- 1
# dummy slurry_prod_rate, wthis should happen inside rates_cpp
slurry_prod_rate <- 2
slurry_prod_rate_vector <- rep(slurry_prod_rate, ncol(cc_fresh))  # or your actual slurry_prod_rate values
names(slurry_prod_rate_vector) <- colnames(cc_fresh)
fresh_derivative <- cc_fresh %*% slurry_prod_rate_vector
fresh_derivative 
####

#### making respiration matrix

cc_resp <- cc_t
cc_resp[,] <- 0

### from earlier calculations on respiration we have this
# note that O2 and H2O are removed because they are not state variables and not tracked
#RFd <- c(RFd = -1, TAN = -0.78, CO2_emis_cum = 2.1, O2 =  -2.1, H2O = 3.44, xa_aer = 0.78,) # 
#starch <- c(starch = -1, TAN = -0.78, CO2_emis_cum = 2.1, O2 =  -2.1, H2O = 3.44, xa_aer = 0.78,) # 
#CPs <- c(CPs = -1, TAN = 0.45725, CO2_emis_cum = 1.28625, O2 = -1.46125, H2O = 0.007249637, xa_aer = 0.54275)
#CPf <- c(CPf = -1, TAN = 0.45725, CO2_emis_cum = 1.28625, O2 = -1.46125, H2O = 0.007249637, xa_aer = 0.54275)
#VSd <-  c(VSd = -1, TAN = -2.24, CO2_emis_cum = 4.7, O2 = -6.685, H2O = 8.69, xa_aer = 2.24)
#Cfat <- c(Cfat = -1, TAN = -9.425, CO2_emis_cum = 3.875, O2 = -25.375, H2O = 39.575, xa_aer = 9.425)

RFd <- c(RFd = -1, TAN = -0.78, CO2_emis_cum = 2.1, xa_aer = 0.78) # 
starch <- c(starch = -1, TAN = -0.78, CO2_emis_cum = 2.1, xa_aer = 0.78) # 
CPs <- c(CPs = -1, TAN = 0.45725, CO2_emis_cum = 1.28625, xa_aer = 0.54275)
CPf <- c(CPf = -1, TAN = 0.45725, CO2_emis_cum = 1.28625, xa_aer = 0.54275)
VSd <-  c(VSd = -1, TAN = -2.24, CO2_emis_cum = 4.7, xa_aer = 2.24)
Cfat <- c(Cfat = -1, TAN = -9.425, CO2_emis_cum = 3.875, xa_aer = 9.425)

cc_resp[names(RFd), 'RFd'] <- RFd
cc_resp[names(starch), 'starch'] <- starch
cc_resp[names(CPs), 'CPs'] <- CPs
cc_resp[names(CPf), 'CPf'] <- CPf
cc_resp[names(VSd), 'VSd'] <- VSd
cc_resp[names(Cfat), 'Cfat'] <- Cfat

# convert to gCOD/mole
moles_to_gCOD <- moles_to_gCOD[colnames(cc_resp)]
cc_resp_COD_mole <- sweep(cc_resp, 1, moles_to_gCOD, `*`)
cc_resp_COD <- cc_resp_COD_mole

cc_resp_COD[,'RFd'] <- -cc_resp_COD_mole[,'RFd']/cc_resp_COD_mole['RFd','RFd']
cc_resp_COD[,'VSd'] <- -cc_resp_COD_mole[,'VSd']/cc_resp_COD_mole['VSd','VSd']
cc_resp_COD[,'starch'] <- -cc_resp_COD_mole[,'starch']/cc_resp_COD_mole['starch','starch']
cc_resp_COD[,'CPs'] <- -cc_resp_COD_mole[,'CPs']/cc_resp_COD_mole['CPs','CPs']
cc_resp_COD[,'CPf'] <- -cc_resp_COD_mole[,'CPf']/cc_resp_COD_mole['CPf','CPf']
cc_resp_COD[,'Cfat'] <- -cc_resp_COD_mole[,'Cfat']/cc_resp_COD_mole['Cfat','Cfat']

# sub_resp in rates_cpp dummy (the state vars that can be respired) 
sub_resp <- sum(y[c('RFd','starch','CPs','CPf','VSd','Cfat')]) # units are gCOD, g , gN, gS
# respiration is also made inside rates_cpp, because depends on temp. Here a dummy
respiration <- 0.4 # unit in rates is gCOD
respiration_sub_resp <- respiration/sub_resp # 

# sweep not neccesary because respiration_sub_resp is the same for all compounds!
# so we can just do multiplication
cc_resp_COD <- cc_resp_COD * respiration_sub_resp
# now make the derivative with the dot product with state vars.
resp_derivative <- cc_resp_COD %*% y[y_names]
resp_derivative
# check the subpart derivatives if they match
print(resp_derivative + alpha_derivative + fresh_derivative)

### NOW the methanogenesis and sulfate reduction! 

predFerm('CH3COOH')
library(biogas)


yield_meth = 0.05
acefrac = 0.66
# can be calculated in ABM, because yield is defined there. 
CH3COOH <- predMethan('CH3COOH', fs = yield_meth)/(predMethan('CH3COOH', fs = yield_meth)[['CH3COOH']]) * acefrac * -1
CH3COOH['H2'] <- 0
H2 <- predMethan('H2', fs = yield_meth)/(predMethan('H2', fs = yield_meth)[['H2']]) * (1-acefrac) * -1
H2['CH3COOH'] <- 0
CH3COOH <- CH3COOH[names(H2)]
CH3COOH

VFA_meth <- CH3COOH + H2
VFA_meth['H2'] <- VFA_meth['H2'] * biogas::molMass('H2') * biogas::calcCOD('H2')
VFA_meth['CO2_emis_cum'] <- VFA_meth['CO2'] * 44.01
VFA_meth['CH3COOH'] <- VFA_meth['CH3COOH'] * biogas::molMass('CH3COOH') * biogas::calcCOD('CH3COOH')
VFA_meth['TAN'] <- VFA_meth['NH3'] * 14.007
VFA_meth['xa'] <- VFA_meth["C5H7O2N"] * biogas::molMass('C5H7O2N') * biogas::calcCOD('C5H7O2N')
VFA_meth['CH4_emis_cum'] <- VFA_meth['CH4'] * biogas::molMass('CH4')
VFA_meth['VFA'] <- VFA_meth['CH3COOH'] + VFA_meth['H2']
# remove names
VFA_meth <- VFA_meth[!names(VFA_meth) %in% c('H2','CH3COOH', 'C5H7O2N','H2O','CO2','NH3','CH4')]
# convert to per gCOD of VFA_meth consumed
VFA_meth <- -VFA_meth/VFA_meth['VFA']

yield_sr = 0.065

CH3COOH_sr_cell <- c(CH3COOH = -1, NH3 = -8/20, C5H7O2N = 8/20, H2O = 1.2) * yield_sr  
CH3COOH_sr_energy <- c(SO4 = -1, CH3COOH = -1, H2O = 2, HS = 1, CO2 = 2) * (1 - yield_sr)

all_names <- union(names(CH3COOH_sr_cell), names(CH3COOH_sr_energy))

CH3COOH_sr_total <- setNames(
  (ifelse(all_names %in% names(CH3COOH_sr_cell), CH3COOH_sr_cell[all_names], 0)) +
    (ifelse(all_names %in% names(CH3COOH_sr_energy), CH3COOH_sr_energy[all_names], 0)),
  all_names
)

CH3COOH_sr_total

H2_sr_cell <- c(H2 = -0.5, NH4 = -1/20, HCO3 = -1/20, C5H7O2N = 1/20, H2O = 9/20) * yield_sr  
H2_sr_energy <- c(SO4 = -1/8, H = -3/16, H2 = -0.5, H2S = 1/16, HS = 1/16, H2O = 0.5) * (1 - yield_sr)

all_names <- union(names(H2_sr_cell), names(H2_sr_energy))

H2_sr_total <- setNames(
  (ifelse(all_names %in% names(H2_sr_cell), H2_sr_cell[all_names], 0)) +
    (ifelse(all_names %in% names(H2_sr_energy), H2_sr_energy[all_names], 0)),
  all_names
)

H2_sr_total <- -H2_sr_total/H2_sr_total['H2']

all_names <- union(names(H2_sr_total), names(CH3COOH_sr_total))

VFA_sr <- setNames(
  (ifelse(all_names %in% names(H2_sr_total), H2_sr_total[all_names], 0)) * (1-acefrac) +
    (ifelse(all_names %in% names(CH3COOH_sr_total), CH3COOH_sr_total[all_names], 0)) * acefrac,
  all_names
)

# change names in VFA_sr to fit to matrix. Look above for methanogens
VFA_sr['H2'] <- VFA_sr['H2'] * biogas::molMass('H2') * biogas::calcCOD('H2')
VFA_sr['CO2_emis_cum'] <- (VFA_sr['CO2'] + VFA_sr['HCO3']) * 44.01
VFA_sr['CH3COOH'] <- VFA_sr['CH3COOH'] * biogas::molMass('CH3COOH') * biogas::calcCOD('CH3COOH')
VFA_sr['TAN'] <- VFA_sr['NH3'] * 14.007
VFA_sr['xa'] <- VFA_sr["C5H7O2N"] * biogas::molMass('C5H7O2N') * biogas::calcCOD('C5H7O2N')
VFA_sr['VFA'] <- VFA_sr['CH3COOH'] + VFA_sr['H2']
VFA_sr['sulfate'] <- VFA_sr['SO4'] * biogas::molMass('S')
VFA_sr['sulfide'] <- (VFA_sr['HS'] + VFA_sr['H2S']) * biogas::molMass('S')
# remove names
VFA_sr <- VFA_sr[!names(VFA_sr) %in% c('H2','CH3COOH', 'C5H7O2N','H2O','CO2','NH3','NH4','CH4','HCO3','HS','H2S','SO4','H')]
# convert to per gCOD of VFA_sr consumed
VFA_sr <- VFA_sr/VFA_sr['VFA']*-1

# check the total stoichiometry of methanogenesis and sulfate reduction
print(VFA_meth)
print(VFA_sr)

# dummy for rut and positions!
rut <- c(m0 = 0.01, m1 = 0.02, m2 = 0.03, sr1 = 0.015)
i_meth <- c(1,2,3)
i_sr <- 4

# edit names for merge consistency
all_names <- union(names(VFA_meth), names(VFA_sr))
VFA_meth[setdiff(all_names, names(VFA_meth))] <- 0
VFA_sr[setdiff(all_names, names(VFA_sr))] <- 0

# apply same order of names
VFA_meth <- VFA_meth[names(VFA_sr)]
VFA <- (VFA_meth * sum(rut[i_meth]) + VFA_sr * sum(rut[i_sr]))

# need to have each xa in this matrix, but for now it is just the sum.
cc_meth_sr <- cc_t
cc_meth_sr[,] <- 0

# add to matrix
cc_meth_sr[names(VFA), 'VFA'] <- VFA

state_vector  <- rep(1, ncol(cc_meth_sr))  # or your actual slurry_prod_rate values
names(state_vector) <- colnames(cc_meth_sr)

rut_derivative <- cc_meth_sr %*% state_vector
print(rut_derivative)

total_derivative <- rut_derivative + fresh_derivative + resp_derivative + alpha_derivative

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
y_names_alpha <- c('xa_dead','xa_bac','RFd','VSd','starch','CPs','CPf','Cfat','VFA','urea','TAN','CO2_emis_cum')

# calc molar stoichiometry for compounds being hydrolyzed through alpha
# xa_dead is tricky because we have tried to separate it from Nitrogen part. All the COD of xa_dead there goes into VFA here.
xa_dead <- c(CH3COOH = 1, xa_dead = -1)

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
names_stoich <- c(names(RFd), 'starch', 'CPs', 'CPf','Cfat','VSd','urea', 'xa_dead')

# make matrix with state vars in rows and stoich in cols
cc <- matrix(nrow = length(y_names_alpha), ncol = length(names_stoich), dimnames = list(y_names_alpha, names_stoich))

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

# convert moles to gCOD
# C5H7O2N = 113.113 g/mol * 1.414515 gCOD/g
moles_to_gCOD <- c(H2O = 16, 
                   starch = biogas::molMass('C6H10O5') * biogas::calcCOD('C6H10O5'), 
                   RFd = biogas::molMass('C6H10O5') * biogas::calcCOD('C6H10O5'), 
                   TAN = 14.007, 
                   H2 = biogas::molMass('H2') * biogas::calcCOD('H2'), 
                   CO2_emis_cum = 44.01, 
                   CH3COOH = biogas::molMass('CH3COOH') * biogas::calcCOD('CH3COOH'), 
                   xa_bac = biogas::molMass('C5H7O2N') * biogas::calcCOD('C5H7O2N'), 
                   xa_dead = biogas::molMass('C5H7O2N') * biogas::calcCOD('C5H7O2N'),
                   CPs = biogas::molMass('C4H6.1O1.2N') * biogas::calcCOD('C4H6.1O1.2N'), 
                   CPf = biogas::molMass('C4H6.1O1.2N') * biogas::calcCOD('C4H6.1O1.2N'), 
                   Cfat = biogas::molMass('C51H98O6') * biogas::calcCOD('C51H98O6'), 
                   VSd = biogas::molMass('C15.9H26.34O9.2N') * biogas::calcCOD('C15.9H26.34O9.2N'),
                   urea = 28.014)

# reorder to ensure same order
moles_to_gCOD <- moles_to_gCOD[colnames(cc)]
# check consistency in lengths
length(moles_to_gCOD) == ncol(cc)
# change unit of matrix to gCOD/mol or gN/mol or g/mol. 
# Need to convert now because we need to merge H2 and CH3COOH and they 
# should be in COD units before doign that. 
cc_COD <- sweep(cc, 2, moles_to_gCOD, `*`)
# remove water and rename relevant column names
cc_COD <- cc_COD[, colnames(cc_COD) != 'H2O']
# H2 + CH3COOH = VFA column
cc_COD <- cbind(cc_COD, VFA = cc_COD[,'H2'] + cc_COD[, 'CH3COOH']) 
# remove H2 and CH3COOH
cc_COD <- cc_COD[, !colnames(cc_COD) %in% c('H2','CH3COOH')]
# reorder column names in same order as rows
cc_COD <- cc_COD[, y_names_alpha]
# transpose it, such that rows = state variable and columns are coefficients multiplied with column state var in dot-product
cc_t <- t(cc_COD)

# cc_t should be passed into rates environment
# alpha should be calculated and depends on e.g. 
# temperature, which is dynamic and has to be in rates
# creating alpha dummy. Als for compounds not being hydrolyzed. It does
# not matter because they have no coefficients anyway (see VFA or bac column)
alpha <- rep(1, length(y_names_alpha)) 
# multiplying with alpha along each column
cc_alpha <- sweep(cc_t, 2, alpha, `*`)
# this is the change in the state variables when only accounting for hydrolysis
# and since y[y_names_alpha] has units of moles here.
# Ideally should hard code the stoich in units of COD instead 
# When 1 gCOD is consumed how much is produced of the stoichiometry 
alpha_derivative <-  cc_alpha %*% y[y_names_alpha] 






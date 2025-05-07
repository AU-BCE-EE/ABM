# Concept for stoich()
# S. Hafner & F. Dalby

# Substrate state variables
source('../../R/predFerm.R')

source('../../R/readFormula.R')

predFerm('C6H12N2O6', acefrac = 0.66, fs = 0.1175)

y <- c(xa = 0,
       slurry_mass = 0, 
       xa_aer = 1,
       xa_bac = 1,
       xa_dead = 1, 
       RFd = 1,
       iNDF = 1,
       ash = 1,
       VSd = 1,
       starch = 1,
       CPs = 1,
       CPf = 1,
       Cfat = 1,
       VFA = 1, 
       urea = 1, 
       TAN = 1, 
       sulfate = 1, 
       sulfide = 1)

RFd <- predFerm('C6H10O5', acefrac = 0.66, fs = 0.1175)
starch <- predFerm('C6H10O5', acefrac = 0.66, fs = 0.1175)
CPs <- predFerm('C4H6.1O1.2N', acefrac = 0.66, fs = 0.1175)
CPf <- predFerm('C4H6.1O1.2N', acefrac = 0.66, fs = 0.1175)
VSd <- predFerm('C15.9H26.34O9.2N', acefrac = 0.66, fs = 0.1175)
Cfat <- predFerm('C51H98O6', acefrac = 0.66, fs = 0.1175)

names_stoich <- c(names(RFd), 'C4H6.1O1.2N','C51H98O6', 'C15.9H26.34O9.2N')

cc <- matrix(nrow = length(y), ncol = length(names_stoich), dimnames = list(names(y), names_stoich))
cc['RFd', names(RFd)] <- RFd
cc['starch', names(starch)] <- starch
cc['CPs', names(CPs)] <- CPs
cc['CPf', names(CPf)] <- CPf
cc['VSd', names(VSd)] <- VSd
cc['Cfat', names(Cfat)] <- Cfat
cc['xa_bac' names(xa)]


# Coef matrix, set in parameters, or calculated before lsoda() call e.g., if composition of VSd varies
cc <- matrix(rep(1, 8), nrow = 4, dimnames = list(names(y), c('CO2', 'VFA')))



# Hydrolysis rates, all 0.1 here
alpha <- y 
alpha[] <- 0.1
cc

# Derivatives 
# Substrate state variables
dydt <- alpha * y # + influent

# Products state variables
dcdt <- colSums(cc * dydt)

dydt
dcdt

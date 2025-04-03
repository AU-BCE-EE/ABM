# Fermentation stoichiometry
# Example calls:
# source('readFormula.R')
# predFerm('C6H10O5', acefrac = 0, fs = 0.1)
# predFerm('C6H10O5', acefrac = 0.5, fs = 0.1)
# predFerm('C6H10O5', acefrac = 1)
# predFerm('C6H10O5', acefrac = 1)

# Function to get stoichiometry for custom organic reaction (O-19)
customOrgStoich <- function(
  form, 
  elements =  c('C', 'H', 'O', 'N')
  ) {
  
  fc <- readFormula(form, elements)

  # Use symbols from O-19 in R&M
  n <- as.numeric(fc['C'])
  a <- as.numeric(fc['H'])
  b <- as.numeric(fc['O'])
  cc <- as.numeric(fc['N'])
  d <- 4 * n + a - 2 * b - 3 * cc
  
  # Put together
  #rr <- c(CO2 = - (n - cc) / d, NH4. = - cc / d, HCO3. = - cc / d, H. = -1, H2O = (2*n - b + cc) /d)
  rr <- c(CO2 = - n / d, NH3 = - cc / d, H. = -1, H2O = (2*n - b + 0*cc) /d)
  rr[form] <-  1/d
  
  return(rr)
  
}

predFerm <- function(
  subform = NULL,           # Character chemical formula of substrate
  biomassform = 'C5H7O2N',  # Biomass empirical formula
  acefrac = 0.5,            # Acetate (vs. H2) fraction
  fs = 0,                    # Fraction substrate going to cell synthesis, fs in Rittmann and McCarty
  elements = c('C', 'H', 'O', 'N'),
  order = 'sort',
  drop = TRUE,
  tol = 1E-10
  ) {

  # Donor half reaction
  rd <- customOrgStoich(subform, elements = elements)

  # Synthesis half reaction
  rc <- customOrgStoich(biomassform, elements = elements)

  # Acceptor reactions
  # Acetate production
  #raa <- c(CO2 = - 1/8, HCO3. = - 1/8, H. = -1, CH3COO. = 1/8, H2O = 3/8)
  raa <- c(CO2 = - 1/4, H. = -1, CH3COOH = 1/8, H2O = 1/4)
  # Hydrogen production
  rah <- c(H. = - 1, H2 = 1 / 2)

  ii <- unique(names(c(rd, rc, raa, rah)))

  # Blanks
  rd[ii[!ii %in% names(rd)]] <- 0
  rc[ii[!ii %in% names(rc)]] <- 0
  raa[ii[!ii %in% names(raa)]] <- 0
  rah[ii[!ii %in% names(rah)]] <- 0

  # Order
  rd <- rd[ii]
  rc <- rc[ii]
  raa <- raa[ii]
  rah <- rah[ii]

  # Acceptor reaction
  ra <- acefrac * raa + (1 - acefrac) * rah 

  fe <- 1 - fs
  
  # Combine
  rtot <- fe * ra + fs * rc  - rd

  ## Simplify: Add TIC elements, combine NH4+ and H+
  #rtot['H.'] <- rtot['H.'] - rtot['HCO3.'] - rtot['CH3COO.'] + rtot['NH4.']
  #rtot['CO2'] <- rtot['CO2'] + rtot['HCO3.']
  #rtot['NH3'] <- rtot['NH4.']
  #rtot['CH3COOH'] <- rtot['CH3COO.'] 
  #rtot['CH3COO.'] <- rtot['NH4.'] <- rtot['HCO3.'] <- 0

  rtot[abs(rtot) < tol] <- 0 
  
  # Drop empty elements
  if (drop) {
    rtot <- rtot[rtot != 0]
  }

  if (!is.na(order[1]) && tolower(order[1]) == 'sort') {
    rtot <- rtot[order(rtot < 0, abs(rtot), decreasing = TRUE)]
  } else if (!is.na(order[1]) && all(sort(order) == sort(names(rtot)))) {
    rtot <- rtot[order]
  } else if (!is.na(order[1])) {
    warning('order argument ignored')
  }

  return(rtot)

}

# For methanogenesis, only difference from ferm is acceptor reaction
# These should be combined--too much code copied!
predMethan <- function(
  subform = NULL,           # Character chemical formula of substrate
  biomassform = 'C5H7O2N',  # Biomass empirical formula
  fs = 0,                   # Fraction substrate going to cell synthesis, fs in Rittmann and McCarty
  elements = c('C', 'H', 'O', 'N'),
  order = 'sort',
  drop = TRUE,
  tol = 1E-10
  ) {

  # Donor half reaction
  rd <- customOrgStoich(subform, elements = elements)

  # Synthesis half reaction
  rc <- customOrgStoich(biomassform, elements = elements)

  # Acceptor reactions
  ra <- c(CO2 = - 1/8, H. = -1, CH4 = 1/8, H2O = 1/4)

  ii <- unique(names(c(rd, rc, ra)))

  # Blanks
  rd[ii[!ii %in% names(rd)]] <- 0
  rc[ii[!ii %in% names(rc)]] <- 0
  ra[ii[!ii %in% names(ra)]] <- 0

  # Order
  rd <- rd[ii]
  rc <- rc[ii]
  ra <- ra[ii]

  fe <- 1 - fs
  
  # Combine
  rtot <- fe * ra + fs * rc  - rd

  ## Simplify: Add TIC elements, combine NH4+ and H+
  #minsp <- c('H.', 'HCO3.', 'CH3COOH', 'CH3COO.', 'NH4.', 'CO2', 'NH3', subform, biomassform)
  #rtot[minsp[!minsp %in% names(rtot)]] <- 0
  #rtot['H.'] <- rtot['H.'] - rtot['HCO3.'] - rtot['CH3COO.'] + rtot['NH4.']
  #rtot['CO2'] <- rtot['CO2'] + rtot['HCO3.']
  #rtot['NH3'] <- rtot['NH4.']
  #rtot['CH3COOH'] <- rtot['CH3COOH'] + rtot['CH3COO.'] 
  #rtot['CH3COO.'] <- rtot['NH4.'] <- rtot['HCO3.'] <- 0

  rtot[abs(rtot) < tol] <- 0 
  
  # Drop empty elements
  if (drop) {
    rtot <- rtot[rtot != 0]
  }

  if (!is.na(order[1]) && tolower(order[1]) == 'sort') {
    rtot <- rtot[order(rtot < 0, abs(rtot), decreasing = TRUE)]
  } else if (!is.na(order[1]) && all(sort(order) == sort(names(rtot)))) {
    rtot <- rtot[order]
  } else if (!is.na(order[1])) {
    warning('order argument ignored')
  }

  return(rtot)

}

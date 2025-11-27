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

# General Rittmann and McCarty stype stoichiometry calculations
RMStoich <- function(subform, rd, ra, rc, fs, dropzero, dropsub, order, tol) {

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

  # Drop substrate if requested
  if (isTRUE(dropsub)) {
    # Adjust coefficients to 1 mol substrate
    rtot <- - rtot / rtot[subform]
    rtot <- rtot[names(rtot) != subform] 
  }

  rtot[abs(rtot) < tol] <- 0 
  
  # Drop empty elements
  if (dropzero) {
    rtot <- rtot[rtot != 0]
  }

  if (!is.na(order[1]) && tolower(order[1]) == 'sort') {
    rtot <- rtot[order(rtot < 0, abs(rtot), decreasing = TRUE)]
  } else if (!is.na(order[1]) && all(sort(order) == sort(names(rtot)))) {
    rtot <- rtot[order]
  } else if (inherits(order, 'logical') && isFALSE(order)) {
    # Skips order
  } else if (!is.na(order[1])) {
    warning('order argument ignored')
  }

  return(rtot)

}

predFerm <- function(
  subform = NULL,           # Character chemical formula of substrate
  biomassform = 'C5H7O2N',  # Biomass empirical formula
  acefrac = 0.5,            # Acetate (vs. H2) fraction
  fs = 0,                    # Fraction substrate going to cell synthesis, fs in Rittmann and McCarty
  elements = c('C', 'H', 'O', 'N'),
  order = 'sort',
  dropzero = TRUE,
  dropsub = FALSE,
  tol = 1E-10
  ) {

  if (biomassform %in% subform) {
    stop('The biomassform formula is the/a substrate! Change one (change in just order is OK).')
  }

  # Vectorize, return matrix without substrate (for 1 mole substrate)
  if (length(subform) > 1) {
    res <- lapply(as.list(subform), predFerm, biomassform = biomassform, acefrac = acefrac, 
                  fs = fs, elements = elements, order = FALSE, dropzero = FALSE, dropsub = TRUE, 
                  tol = tol)
    resmat <- matrix(unlist(res), nrow = length(subform), byrow = TRUE)
    colnames(resmat) <- names(res[[1]])
    rownames(resmat) <- subform

    if (dropzero) {
      resmat <- resmat[, colSums(resmat) != 0]
    }
    return(resmat)
  }

  # Donor half reaction
  rd <- customOrgStoich(subform, elements = elements)
  
  ## Donor needs H2 and CH3COOH
  #for (sp in c('H2', 'CH3COOH')) {
  #  if (! sp %in% names(rd)) {
  #    rd[sp] <- 0
  #  }
  #}
  #rd <- rd[sort(names(rd))]

  # Synthesis half reaction
  rc <- customOrgStoich(biomassform, elements = elements)

  # Acceptor reactions
  # Acetate production
  #raa <- c(CO2 = - 1/8, HCO3. = - 1/8, H. = -1, CH3COO. = 1/8, H2O = 3/8)
  raa <- c(CO2 = - 1/4, H. = -1,  CH3COOH = 1/8, H2O = 1/4, H2 = 0)
  # Hydrogen|production |            |            |          |
  rah <- c(CO2 = 0,     H. = - 1, CH3COOH = 0,   H2O = 0,   H2 = 1/2)
  
  # Acceptor reaction
  ra <- acefrac * raa + (1 - acefrac) * rah 

  rtot <- RMStoich(subform = subform, rd = rd, ra = ra, rc = rc, fs = fs, 
                   dropzero = dropzero, dropsub = dropsub, order = order, 
                   tol = tol)

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
  dropzero = TRUE,
  dropsub = FALSE,
  tol = 1E-10
  ) {

  if (biomassform %in% subform) {
    stop('The biomassform formula is the/a substrate! Change one (change in just order is OK).')
  }

  # Vectorize, return matrix without substrate (for 1 mole substrate)
  if (length(subform) > 1) {
    res <- lapply(as.list(subform), predMethan, biomassform = biomassform, 
                  fs = fs, elements = elements, order = FALSE, dropzero = FALSE, dropsub = TRUE, 
                  tol = tol)
    resmat <- matrix(unlist(res), nrow = length(subform), byrow = TRUE)
    colnames(resmat) <- names(res[[1]])
    rownames(resmat) <- subform

    if (dropzero) {
      resmat <- resmat[, colSums(resmat) != 0]
    }
    return(resmat)
  }

  # Donor half reaction
  rd <- customOrgStoich(subform, elements = elements)

  ## Donor needs H2 and CH3COOH
  #for (sp in c('H2', 'CH3COOH')) {
  #  if (! sp %in% names(rd)) {
  #    rd[sp] <- 0
  #  }
  #}
  #rd <- rd[sort(names(rd))]

  # Synthesis half reaction
  rc <- customOrgStoich(biomassform, elements = elements)

  # Acceptor reaction
  ra <- c(CO2 = - 1/8, H. = -1, CH4 = 1/8, H2O = 1/4)

  rtot <- RMStoich(subform = subform, rd = rd, ra = ra, rc = rc, fs = fs, 
                   dropzero = dropzero, dropsub = dropsub, order = order, 
                   tol = tol)

  return(rtot)

}

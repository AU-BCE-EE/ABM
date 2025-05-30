stoich_species2 <- function(acefrac, fs, y) {
  
  formulas <- c(starch = 'C6H10O5', 
                Cfat = 'C51H98O6', 
                CPs = 'C4H6.1O1.2N',
                CPf = 'C4H6.1O1.2N',
                RFd = 'C6H10O5',
                VSd = 'C15.9H26.34O9.2N')
  
  ferm_results <- lapply(formulas, function(f) predFerm(f, acefrac = acefrac, fs = fs))
  ferm_ref <- sapply(seq_along(formulas), function(i) ferm_results[[i]][[formulas[i]]])
  
  # normalize to COD units
  stoich_list <- mapply(function(res, ref) res / ref * -1, ferm_results, ferm_ref, SIMPLIFY = FALSE)
  
  # rename primary compound
  for (nm in names(stoich_list)) {
    names(stoich_list[[nm]])[names(stoich_list[[nm]]) == formulas[nm]] <- nm
  }
  
  # fixed stoichiometries
  xa_dead <- c(xa_dead = -1, CH3COOH = 1)
  urea <- c(urea = -1, NH3 = 1, CO2 = 0.5, H2O = -0.5)
  
  # complete stoichiometry names
  names_stoich <- c(names(stoich_list$RFd), 'starch', 'CPs', 'CPf','Cfat', 'VFA', 'VSd','urea', 'xa_dead','sulfide','sulfate','xa','slurry_mass',
                    'xa_aer','iNDF','ash','NH3_emis_cum','N2O_emis_cum','CH4_emis_cum','COD_conv_cum','COD_conv_cum_meth',
                    'COD_conv_cum_respir','COD_conv_cum_sr', 'COD_load_cum','C_load_cum', 'N_load_cum', 'slurry_load_cum')
  
  # create empty matrix
  cc <- matrix(0, nrow = length(y), ncol = length(names_stoich), dimnames = list(names(y), names_stoich))
  
  # fill matrix with stoichiometries
  for (nm in names(stoich_list)) cc[nm, names(stoich_list[[nm]])] <- stoich_list[[nm]]
  cc['xa_dead', names(xa_dead)] <- xa_dead
  cc['urea', names(urea)] <- urea
  
  # rename columns
  colnames(cc)[colnames(cc) == "NH3"] <- "TAN"
  colnames(cc)[colnames(cc) == "C5H7O2N"] <- "xa_bac"
  colnames(cc)[colnames(cc) == "CO2"] <- "CO2_emis_cum"
  
  # precompute molar masses and COD equivalents only once
  compounds <- c("H2O","C6H10O5","C51H98O6","C4H6.1O1.2N","C15.9H26.34O9.2N",
                 "H2","CH3COOH","C5H7O2N")

  mol_masses <- sapply(compounds, biogas::molMass)
  cod_vals   <- sapply(compounds, biogas::calcCOD)
  
  moles_to_gCOD <- c(
    H2O = 16,
    starch = mol_masses["C6H10O5"] * cod_vals["C6H10O5"],
    RFd = mol_masses["C6H10O5"] * cod_vals["C6H10O5"],
    TAN = 14.007,
    H2 = mol_masses["H2"] * cod_vals["H2"],
    CO2_emis_cum = 44.01,
    CH3COOH = mol_masses["CH3COOH"] * cod_vals["CH3COOH"],
    xa_bac = mol_masses["C5H7O2N"] * cod_vals["C5H7O2N"],
    xa_aer = mol_masses["C5H7O2N"] * cod_vals["C5H7O2N"],
    xa_dead = mol_masses["C5H7O2N"] * cod_vals["C5H7O2N"],
    xa = mol_masses["C5H7O2N"] * cod_vals["C5H7O2N"],
    CPs = mol_masses["C4H6.1O1.2N"] * cod_vals["C4H6.1O1.2N"],
    CPf = mol_masses["C4H6.1O1.2N"] * cod_vals["C4H6.1O1.2N"],
    Cfat = mol_masses["C51H98O6"] * cod_vals["C51H98O6"],
    VSd = mol_masses["C15.9H26.34O9.2N"] * cod_vals["C15.9H26.34O9.2N"],
    urea = 28.014
  )
  
  moles_to_gCOD_clean <- rep(0, length(colnames(cc)))
  names(moles_to_gCOD_clean) <- colnames(cc)
  
  # Now, fill in the known values
  known_values <- c(
    H2O = 16, 
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
    urea = 28.014
  )
  
  # Populate the matching entries
  moles_to_gCOD_clean[names(known_values)] <- known_values
  
  # convert coefficients to gCOD/mol
  cc_gCOD_mole <- sweep(cc, 2, moles_to_gCOD_clean, `*`)
  
  # normalize each degradation compound row by itself
  for (comp in c("xa_dead","RFd","VSd","starch","CPs","CPf","Cfat","urea")) {
    cc_gCOD_mole[comp, ] <- -cc_gCOD_mole[comp, ] / cc_gCOD_mole[comp, comp]
  }
  
  # H2 + CH3COOH → VFA
  cc_gCOD_mole[, 'VFA'] <- cc_gCOD_mole[,'H2'] + cc_gCOD_mole[, 'CH3COOH']
  
  # remove H2, CH3COOH, H2O
  drop_cols <- intersect(c("H2","CH3COOH","H2O"), colnames(cc_gCOD_mole))
  cc_gCOD_mole <- cc_gCOD_mole[, !colnames(cc_gCOD_mole) %in% drop_cols]
  
  # reorder columns to match y
  cc_gCOD_mole <- cc_gCOD_mole[, names(y)]
  
  return(t(cc_gCOD_mole))
}
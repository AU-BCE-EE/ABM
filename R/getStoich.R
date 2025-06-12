# Figure out stoichiometry matrix of fermentation substrates from elemental formula

getStoich <- function(pars) {

  # Get molar stoichiometric coefficients
  st <- predFerm(pars$forms)

  # Drop water (ignored, treated as conservative)
  st <- st[rownames(st) != 'H2O', ] 

  # Switch from chemical formulas of columns to substrate names
  colnames(st) <- names(pars$forms)

  # Switch to master species names for products (some match)
  rownames(st) <- pars$mspec[rownames(st)]

  # Adjust coefficients to COD mass, N mass, C mass, S mass, or total mass
  for (i in 1:nrow(st)) {
    ff <- rownames(st)[i]
    st[i, ] <- st[i, ] * 1 / pars$mcf[ff]
  }
  
  for (i in 1:ncol(st)) {
    ff <- colnames(st)[i]
    st[, i] <- st[, i] * pars$mcf[ff]
  }

  return(st)

}

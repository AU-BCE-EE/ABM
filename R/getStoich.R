# Figure out stoichiometry matrix of fermentation substrates from elemental formula

getStoich <- function(pars) {

  # Get molar stoichiometric coefficients
  # Vectorize, return matrix without substrate (for 1 mole substrate)
  res <- lapply(as.list(pars$forms), predFerm)
  # Align names and sort before combining in matrix
  nn <- unique(unlist(lapply(res, names)))
  for (i in 1:length(res)) {
    res[[i]][nn[!nn %in% names(res[[i]])]] <- 0
    res[[i]] <- res[[i]][nn]
  }
  st <- matrix(unlist(res), ncol = length(pars$forms), byrow = FALSE)
  rownames(st) <- names(res[[1]])
  colnames(st) <- pars$forms

  # Drop 0
  st <- st[rowSums(st) != 0, , drop = FALSE]
  
  # Drop water (ignored, treated as conservative)
  st <- st[rownames(st) != 'H2O', , drop = FALSE] 

  # Switch from chemical formulas of columns to substrate names
  colnames(st) <- names(pars$forms)

  # Switch to master species names for products (some match)
  rownames(st) <- pars$mspec[rownames(st)]

  # Adjust coefficients to COD mass, N mass, C mass, S mass, or total mass
  for (i in 1:nrow(st)) {
    ff <- rownames(st)[i]
    st[i, ] <- st[i, ] * pars$mcf[ff]
  }
  
  for (i in 1:ncol(st)) {
    ff <- colnames(st)[i]
    st[, i] <- st[, i] * 1 / pars$mcf[ff]
  }

  return(st)

}

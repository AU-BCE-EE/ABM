# Function for emptying storage struction
# Not exported, only called by abm_*

# Note: enrich_names could also be length 1 vector: '^xa|^RFd$|^iNDF$' etc.

emptyStore <- function(y, resid_mass, resid_enrich, 
                       enrich_names = c('^xa', '^m[0-9]', '^sr[0-9]', '^RFd$', '^iNDF$', '^VSd$', '^CP$', '^CF$', '^starch$', '^ash$'),
                       ignore_names = c('_cum_', '_conv_', 'cum$')) {

  y <- unlist(y)
  slurry_mass <- y['slurry_mass']

  which.ignore <- grepl(paste(ignore_names, collapse = '|'), names(y))

  if (slurry_mass > resid_mass) {
    # Masses before emptying
    y.before <- y
    # Calculate mass of each variable remaining in storage
    resid_frac <- resid_mass / slurry_mass
    resid_par <- logistic(logit(resid_frac) + resid_enrich)
    which.enrich <- grepl(paste(enrich_names, collapse = '|'), names(y))
    y[which.enrich] <- y[which.enrich] * resid_par
    y[!which.enrich & !which.ignore] <- y[!which.enrich & !which.ignore] * resid_frac
    # Effluent
    y.eff <- y.before - y
    y.eff <- y.eff[!which.ignore]
    names(y.eff) <- paste0(names(y.eff), '_eff')
    y <- c(y, y.eff)
  } else {
    warning('Emptying skipped because of low slurry level.')
    y.eff <- 0 * y
    y.eff <- y.eff[!which.ignore]
    names(y.eff) <- paste0(names(y.eff), '_eff')
    y <- c(y, y.eff)
  }

  return(y)

}

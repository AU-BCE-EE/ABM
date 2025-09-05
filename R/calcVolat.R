
calcVolat <- function(p, vi) {

  volat <- 0 * vi

  # Emission of other species                                           
  if (!is.null(p$kl)) {
    # Check for missing species (should probably move, but this is better than useless error)
    if (any(is.na(p$conc_sp[names(p$kl)]))) {
      cat('Seems like a species definition is missing. . .')
      browser()
      stop('Seems like a species definition is missing. . .')
    }
    volat[paste0(names(p$kl), '_emis_cum')] <- p$kl * p$conc_sp[names(p$kl)] * p$area
    # Remove emitted amount from component pool
    volat[p$mspec[names(p$kl)]] <- - volat[paste0(names(p$kl), '_emis_cum')]
  } 

  return(volat)
  
}


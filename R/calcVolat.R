
calcVolat <- function(p, vi) {

  volat <- 0 * vi

  # Emission of other species                                           
  if (!is.null(p$kla)) {
    volat[paste0(names(p$kla), '_emis_cum')] <- p$kla * p$conc_sp[names(p$kla)] * p$area
    # Remove emitted amount from component pool
    volat[p$mspec[names(p$kla)]] <- - volat[paste0(names(p$kla), '_emis_cum')]
  } 

  return(volat)
  
}


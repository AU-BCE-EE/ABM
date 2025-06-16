
calcUt <- function(p, y, qhat, rut) {

  # Sulfate reducers (VFA part and sulfate part)
  if ('SO4m2' %in% names(y)) {
    srvp <- y['CH3COOH'] / (p$ksv[p$srs] * y['slurry_mass'] + y['CH3COOH'])
    srsp <- y['SO4m2'] / (p$kss[p$srs] * y['slurry_mass'] + y['SO4m2'])
    srt <- srvp * srsp
    names(srt) <- p$srs
  } else {
    srt <- (rep(0, length(p$srs)))
    names(srt) <- p$srs
  }
  
  # Methanogens (VFA part is the only part)
  met <- y['CH3COOH'] / (p$ksv[p$meths] * y['slurry_mass'] + y['CH3COOH'])

  # Substrate utilization rate
  rut[p$grps] <- p$ired[p$grps] * qhat[p$grps] * y[p$grps] * c(met, srt)[p$grps]

  if ('SO4m2' %in% names(y)) {
    rut['SO4m2'] <- - sum(rut[p$srs] - p$yield[p$srs] * rut[p$srs])
    if ('H2S' %in% names(y)) {
      rut['H2S'] <-   sum(srt[p$srs] - p$yield[p$srs] * rut[p$srs])
    }
  }

  return(rut)

}

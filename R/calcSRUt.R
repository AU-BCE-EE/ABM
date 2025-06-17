
calcSRUt <- function(p, y, qhat, rut) {

  # Sulfate reducers (VFA part and sulfate part)
  if (p$sromit) {
    return(rut)
  } 
  
  rut[p$srs] <- p$ired[p$srs] * qhat[p$srs] * y[p$srs] * 
                y['CH3COOH'] / (p$ksv[p$srs] * y['slurry_mass'] + y['CH3COOH']) *
                y['SO4m2'] / (p$kss[p$srs] * y['slurry_mass'] + y['SO4m2'])

  return(rut)

}

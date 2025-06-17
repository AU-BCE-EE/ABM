
calcMethUt <- function(p, y, qhat, rut) {

  # Methanogens (VFA part is the only part)
  halfrut <- y['CH3COOH'] / (p$ksv[p$meths] * y['slurry_mass'] + y['CH3COOH'])
  # Substrate utilization rate
  rut[p$meths] <- p$ired[p$meths] * qhat[p$meths] * y[p$meths] * halfrut[p$meth]

  return(rut)

}

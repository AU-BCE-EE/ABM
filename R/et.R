# Evaporation or evapotranspiration model

et <- function(cmakk = 0.7, temp_C, pres_kpa, rs) {

  # cmakk = C constant in model (0.7 from Hansen 1984)
  # temp_C = average air temperature in degrees C
  # pres_kpa = average air pressure in kPa
  # rs = solar radiation (MJ/m2-d)
  # delta = Slope of vapor pressure vs. temperature (KPa/K)
  # lambda = Latent heat of vaporization (MJ/kg)

  # et = Reference evapotranspiration (mm/d)

  # Set constants
  lambda <- 2.45

  # Get values for constants
  delta <- slopel(temp_C = temp_C)
  gamma <- psycg(pres_kpa = pres_kpa)

  et <- cmakk * delta / (delta + gamma) * rs / lambda 

  return(et)

}

psycg <- function(pres_kpa) {

  # Psychrometric constant kPa/K
  g <- 0.00665 * pres_kpa

  ## Alternative
  #lambda <- 2.45  # Latent heat of vaporization MJ/kg
  #cp <- 0.001013  # Specific heat MJ/kg-K
  #epsilon <- 0.62 # Molecular weight ratio

  #g <- cp * pres_kpa / (epsilon * lambda)

  return(g)
}

slopel <- function(temp_C) {
  # Slope of vapor pressure vs. temperature KPa/K
  d <- 4098 * (0.6108 * exp(17.27 * temp_C / (temp_C + 237.3))) / (temp_C + 237.3)^2

  return(d)
}

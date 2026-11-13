rm(list = ls())

library(data.table)

pH <- 7
temp_K <- 273+35
g_NH4 <- 0.7
slurry_mass <- 1

log_ka <- c(NH3 = - 0.09046 - 2729.31/temp_K, 
            H2S = - 3448.7/temp_K + 47.479 - 7.5227* log(temp_K),
            HAC = -4.8288 + 21.42/temp_K)

NH3_frac <- ((1/(1 + 10^(- log_ka[['NH3']] + log10(g_NH4) - pH))))
NH4_frac <- 1- NH3_frac
HAC_frac <- 1-(1/(1 + 10^(-log_ka[['HAC']] - pH)))
H2S_frac <- 1 - (1/(1 + 10^(- log_ka[['H2S']] - pH))) # H2S fraction of total sulfid

ki_H2S_min <- 0.08
ki_H2S_slope <- -0.1093
ki_H2S_int <- 1.02

if(pH >= 6.8){
  IC50 <- ki_H2S_slope * pH + ki_H2S_int
} else {
  IC50 <- ki_H2S_slope * 6.8 + ki_H2S_int
}

ki_NH3_min <- 0.015
ki_NH3_max <- 0.13
ki_NH4_min <- 2.7
ki_NH4_max <- 4.8

ki_HAC_half <- 0.31

ki_H2S_min = 0.08
ki_H2S_max = 0.4298
ki_HAC_min = 0.05
ki_HAC_half = 0.2625#0.81

NH3_inhib <- NULL
NH4_inhib <- NULL
TAN_conc <- NULL

for(TAN in c(seq(from = 0, to = 20, by = 1/600))){
  
  NH3_inhib1 <- ifelse(NH3_frac * TAN/slurry_mass <= ki_NH3_min, 1, exp(-2.77259 * ((NH3_frac * (TAN/(slurry_mass)) - ki_NH3_min)/(ki_NH3_max - ki_NH3_min))^2))
  NH4_inhib1 <- ifelse((1 - NH3_frac) * (TAN/slurry_mass) <= ki_NH4_min, 1, exp(-2.77259*(((1 - NH3_frac) * (TAN/(slurry_mass)) - ki_NH4_min)/(ki_NH4_max - ki_NH4_min))^2))
  
  NH3_inhib <- c(NH3_inhib, NH3_inhib1)
  NH4_inhib <- c(NH4_inhib, NH4_inhib1)
  TAN_conc <- c(TAN_conc, TAN)
  
}

H2S_inhib_new <- NULL
H2S_inhib <- NULL
H2S_conc <- NULL

for(sulfide in c(seq(from = 0.00, to = 1, by = 1/600))){
  
  a <- -0.5/(IC50-(H2S_frac * ki_H2S_min))
  x <- H2S_frac * sulfide/(slurry_mass)
  b <- 1 -(-0.5/(IC50-(H2S_frac * ki_H2S_min)) * H2S_frac * ki_H2S_min/(slurry_mass))
  ifelse (sulfide/(slurry_mass) > ki_H2S_min, 
          H2S_inhib_meth <- a * x + b, 
          H2S_inhib_meth <- 1)
  
  if(H2S_inhib_meth < 0) {H2S_inhib_meth <- 0}
  
  H2S_inhib <- c(H2S_inhib, H2S_inhib_meth)
  H2S_conc <- c(H2S_conc, sulfide * H2S_frac)
  H2S_inhib_new1 <- ifelse(H2S_frac * sulfide/slurry_mass <= ki_H2S_min, 1, exp(-2.77259 * ((H2S_frac * (sulfide/(slurry_mass)) - ki_H2S_min)/(ki_H2S_max - ki_H2S_min))^2))
  H2S_inhib_new <- c(H2S_inhib_new, H2S_inhib_new1)
  
  
}


HAC_inhib <- NULL
HAC_inhib_new <- NULL
HAC_conc <- NULL

for(VFA in c(seq(from = 0.00, to = 20, by = 1))){
  
  VFA_inhib <- ifelse(HAC_frac * VFA/slurry_mass >= 0.05, 1.16 * ki_HAC_half/(ki_HAC_half + (HAC_frac * VFA/slurry_mass)), 1) 
  
  HAC_inhib <- c(HAC_inhib, VFA_inhib)
  HAC_conc <- c(HAC_conc, VFA * HAC_frac)

  HAC_inhib_new1 <- ifelse(
    HAC_conc >= ki_HAC_min,
    (2 - ki_HAC_half / (ki_HAC_half + ki_HAC_min)) *
      ki_HAC_half / (ki_HAC_half + HAC_conc),
    1.0
  )
  HAC_inhib_new <- c(HAC_inhib_new, HAC_inhib_new1)
  
  
}

HAC_dat_old <- data.frame(conc = HAC_conc, inhib = HAC_inhib, type = "old", comp = "HAC")
H2S_dat_old <- data.frame(conc = H2S_conc, inhib = H2S_inhib, type = "old", comp = "H2S")
HAC_dat_new <- data.frame(conc = HAC_conc, inhib = HAC_inhib_new, type = "new", comp = "HAC")
H2S_dat_new <- data.frame(conc = H2S_conc, inhib = H2S_inhib_new, type = "new", comp = "H2S")

dat <- rbind(HAC_dat_old, H2S_dat_old, HAC_dat_new, H2S_dat_new)

library(ggplot2)

ggplot(HAC_dat_new, aes(conc, inhib, col = type)) + geom_point() + facet_grid(~comp)


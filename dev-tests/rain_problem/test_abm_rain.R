rm(list = ls())
library(data.table)

devtools::load_all()

slurry_mass_dat <- fread('slurry_digestate_dat.csv')


plot(slurry_mass_dat$time, slurry_mass_dat$slurry_mass)

area = 1017.876
rain = 1.9
evap = 0

slurry_mass_dat <- as.data.frame(slurry_mass_dat[, .(time, slurry_mass)])
test <- abm(days = max(slurry_mass_dat$time), add_pars = list(evap = evap, slurry_mass = slurry_mass_dat, rain = rain, area = area))

plot(slurry_mass_dat$time, slurry_mass_dat$slurry_mass)
lines(test$time, test$slurry_mass, col = 'red')
# add some slurry mass to avoid negative slurry mass in start

slurry_mass_dat$slurry_mass <- slurry_mass_dat$slurry_mass + area * rain
test <- abm(days = max(slurry_mass_dat$time), add_pars = list(evap = evap, slurry_mass = slurry_mass_dat, rain = rain, area = area))

rm(list = ls())
#library(ABM)

devtools::load_all('~/GitHub_repos/ABM')

library(data.table)
slurry_mass_dat <- data.frame(fread('slurry_mass.csv'))

# Problem one is with default rain and evap
# First should throw error
#out_a0 <- abm(days = 4*365, add_pars = list(slurry_mass = slurry_mass_dat))
# So adjust area to get reasonable adjusted level (unsure of actual area)
out_a1 <- abm(days = 4*365, add_pars = list(slurry_mass = slurry_mass_dat, area = 100))
# Does it also work with different emptying alignment?
out_a2 <- abm(days = 4*365, add_pars = list(slurry_mass = slurry_mass_dat, area = 100), approx_method = c(temp = 'linear', pH = 'linear', slurry_mass = 'mid')) 
out_a3 <- abm(days = 4*365, add_pars = list(slurry_mass = slurry_mass_dat, area = 100), approx_method = c(temp = 'linear', pH = 'linear', slurry_mass = 'late')) 

out_b <- abm(days = 4*365, add_pars = list(slurry_mass = slurry_mass_dat, rain = 0, evap = 0))

# Check for mismatch between input slurry_mass and output slurry_mass
plot(slurry_mass_dat$time, slurry_mass_dat$slurry_mass)
lines(out_a1$time, out_a1$slurry_mass, col = 'red')
lines(out_a2$time, out_a2$slurry_mass, col = 'orange', lty = 2)
lines(out_a3$time, out_a3$slurry_mass, col = 'purple', lty = 3)
lines(out_b$time, out_b$slurry_mass, col = 'blue')




plot(slurry_mass_dat$time, slurry_mass_dat$slurry_mass, ylim = c(-1E6, 1E6))

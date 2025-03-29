# Sasha's examples for abm() testing during package development


devtools::load_all()

slurry_mass_dat <- data.frame(time = c(0, 25, 50, 100), slurry_mass = c(100, 50, 100, 200))

out1 <- abm(365, add_pars = list(slurry_mass = slurry_mass_dat))

library(ABM)
library(ggplot2)

devtools::load_all('../')

H2SO4_titrat(2, 7, 'pig')


anim_pars <- pig_pars2.0
out <- abm(anim_pars = anim_pars, add_pars = list(slurry_prod_rate = 0, slurry_mass = 5000, evap = 0, rain = 0, temp_C = 20, empty_int = 0, area = 2, floor_area = 0))
out$treatment <- "control"

# g kg/ H2SO4/m3 slurry --> 5 kg S/ m3 slurry
H2SO4_dose <- 3 # kg/m3
anim_pars$conc_fresh$sulfate <- H2SO4_dose * 32/96 * 0.96
anim_pars$pH <- 'calc'

out1 <- abm(anim_pars = anim_pars, add_pars = list(slurry_prod_rate = 0, slurry_mass = 5000, evap = 0, rain = 0, temp_C = 20, empty_int = 0, area = 2, floor_area = 0))
out1$treatment <- "acid"


dat <- rbind(out, out1)

ggplot(dat, aes(time, CH4_emis_rate, col = treatment)) + geom_point() + labs(x = 'days', y = 'CH4, g/day') + theme_minimal()

png('Emis.png', width = 10/2.54, height = 10/2.54, units = 'in', res = 600)
ggplot(dat, aes(time, CH4_emis_rate, col = treatment)) + geom_point() + labs(x = 'days', y = 'CH4, g/day') + theme_minimal()
dev.off()

png('SO4.png', width = 10/2.54, height = 10/2.54, units = 'in', res = 600)
ggplot(dat, aes(time, sulfate_conc, col = treatment)) + geom_point() + labs(x = 'days', y = 'sulfate-S, g/kg') + theme_minimal()
dev.off()

png('pH.png', width = 10/2.54, height = 10/2.54, units = 'in', res = 600)
ggplot(dat, aes(time, pH, col = treatment)) + geom_point() + theme_minimal()
dev.off()




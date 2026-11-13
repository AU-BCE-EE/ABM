
devtools::load_all("..")
library(ggplot2)

out <- cbind(abm(), treatment = 'contr')
outa <- cbind(abm(add_pars = list(pH = 6)), treatment = 'pH 6')

dat <- rbind(out,outa)

ggplot(dat, aes(time, CH4_emis_rate, col = treatment)) + geom_line() +
  geom_point() +
  labs(x = 'Time (days)',
       y = 'CH4 emissions (g/d)',
       col = 'Treatment') +
  theme_minimal()

library(data.table)

setDT(dat)

outb <- cbind(abm(add_pars = list(pH = 6), pH_inhib_overrule = TRUE), treatment = 'pH 6, pH overrule')
dat <- rbind(dat, outb)

ggplot(dat[treatment != 'contr'], aes(time, CH4_emis_rate, col = treatment)) + geom_point() +
  labs(x = 'Time (days)',
       y = 'CH4 emissions (g/d)',
       col = 'Treatment') +
  theme_minimal()

conc_SO4 <- H2SO4_titrat(pH = 6, class_anim = 'pig')
print(paste('the concentration of SO4-S in the fresh manure is approximately', conc_SO4$cSO4))

outc <- cbind(abm(add_pars = list(conc_fresh.sulfate = conc_SO4$cSO4)), treatment = 'SO4 adjust only')
outd <- cbind(abm(add_pars = list(conc_fresh.sulfate = conc_SO4$cSO4, pH = 6)), treatment = 'pH 6 + SO4 adjust')

dat <- rbind(dat, outc, outd)

ggplot(dat, aes(time, CH4_emis_rate, col = treatment)) + geom_point() + facet_grid(~treatment)
labs(x = 'Time (days)',
     y = 'CH4 emissions (g/d)',
     col = 'Treatment') +
  theme_minimal()

dat <- rbind(dat, outc, outd)

dat$rut_meth <- (dat$rut_m0 + dat$rut_m1 + dat$rut_m2)/dat$slurry_mass
dat$rut_sr <- dat$rut_sr1/dat$slurry_mass
dat$VFA_conc_10 <- dat$VFA_conc/10
dat_long <- melt(dat, measure.vars = c('rut_meth','rut_sr','VFA_conc_10'))
ggplot(dat_long, aes(time, value, col = variable)) + geom_point() + facet_grid(~treatment) + 
  labs(x = 'Time (days)',
       y = 'gCOD/kg/day or gCOD/10/kg',
       col = 'type') +
  theme_minimal()

# get the effluent flows from the barn
effluents <- setDT(abm(value = 'eff'))
# define names that we want to extract from the effluents
names_conc <- paste0(names(ABM::man_pars2.0$conc_fresh),'_eff_conc')
# add time to the names
names_conc <- c('time', names_conc)
# get the effluent concentrations which is 
# the fresh concentration for the storage
conc_fresh_dat <- data.frame(effluents[, ..names_conc])
# change names back
names(conc_fresh_dat) <- c('time', names(ABM::man_pars2.0$conc_fresh))
# get slurry mass effluents, which is
# the slurry mass transported to the storage
slurry_mass_dat <- data.frame(effluents[, .(time, slurry_mass = slurry_mass_eff_cum)])

# get some manure storage temperature data
temp_dat <- ABM::outside_slurry_temp_vechi

# run storage simulation without acidification 
oute <- setDT(cbind(abm(add_pars = list(slurry_mass = slurry_mass_dat, 
                                        temp_C = temp_dat, conc_fresh = conc_fresh_dat)), 
                    treatment = 'control storage'))

# what day of year is acid added? Because our temperature data is with reference
# to the Jan 1st our simulation will also be with that reference
# We choose begin = 'January' in the built in doy() function

t_acid <- doy("May", begin = "January")$doy

# We need to calculate the concentration of sulfate that would need to be in 
# the fresh slurry at the time of acidification to reach a pH of 6 in the 
# bulk slurry ion occurs.

# first know slurry mass and slurry production rate at around that time.
slurry_mass_t_acid <- oute[time == t_acid, 'slurry_mass']
slurry_prod_rate <- oute[time == t_acid+1, 'slurry_mass'] - slurry_mass_t_acid 

# calculate the theoretical concentration of sulfate-S in the fresh manure to
# reach the sulfate concentration that gives pH 6 
# Use previous estimate from H2SO4_titrate(), that is conc_SO4$cSO4.
conc_SO4_t_acid <- conc_SO4$cSO4 * slurry_mass_t_acid/slurry_prod_rate

# add to conc_fresh this extra data to the conc_fresh_dat
conc_fresh_dat_t_acid <- conc_fresh_dat[1:3, ]
conc_fresh_dat_t_acid$time <- c(t_acid, t_acid + 1, t_acid + 1.001)
conc_fresh_dat <- rbind(conc_fresh_dat, conc_fresh_dat_t_acid)
conc_fresh_dat <- conc_fresh_dat[order(conc_fresh_dat$time), ]
conc_fresh_dat[conc_fresh_dat$time == t_acid | conc_fresh_dat$time == t_acid+1, 'sulfate'] <- conc_SO4_t_acid

outg <- setDT(cbind(abm(add_pars = list(slurry_mass = slurry_mass_dat, 
                                        temp_C = temp_dat, 
                                        conc_fresh = conc_fresh_dat, pH = 'calc')), 
                    treatment = 'Single dose, 3.6 kg/m3 H2SO4'))
# Lets go further down in pH
conc_SO4_pH_5 <- H2SO4_titrat(pH = 5, class_anim = 'pig')$cSO4
conc_SO4_pH_5 <- conc_SO4_pH_5 * slurry_mass_t_acid/slurry_prod_rate
conc_fresh_dat[conc_fresh_dat$time == t_acid | conc_fresh_dat$time == t_acid+1, 'sulfate'] <- conc_SO4_pH_5

outh <- setDT(cbind(abm(add_pars = list(slurry_mass = slurry_mass_dat, 
                                        temp_C = temp_dat, 
                                        conc_fresh = conc_fresh_dat, pH = 'calc')), 
                    treatment = 'Single dose, 8 kg/m3 H2SO4'))

storage <- rbind(oute, outg, outh)

storage_long <- melt(storage, measure.vars = c('CH4_emis_rate','sulfate_conc','pH'))
# plot the sulfate concentration in the manure
ggplot(storage_long[variable == 'CH4_emis_rate'], aes(time, value, col = treatment)) + 
  geom_point() + 
  facet_wrap(~variable, scales = 'free') + 
  labs(x = 'Time (days)',
       y = 'g CH4/d',
       col = 'Treatment')

plot(outh$time, outh$pH)
plot(outh$time, outh$sr1_conc)


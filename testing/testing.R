
rm(list = ls())
ff <- list.files('../R', full.names = TRUE)
for (i in ff) source(i)

# Define influent composition
# Set concentrations in g/kg
conc_fresh <- list(S2 = 0, SO4 = 0, TAN = 1.0, 
                   TS = 30, TSS = 20, 
                   VS = 20, VSS = 15, 
                   dsVS = 0, dVSS = 10)

# pH of 7 and density of 1000 kg/m3
man_pars <- list(conc_fresh = conc_fresh,
                 pH = 7, dens = 1000)

# Set some other management parameters
# Production (influent) rate in kg/d
mng_pars <- list(slurry_prod_rate = 2E6, slurry_rem_rate = 1E6, 
                 slurry_mass = 100, 
                 storage_depth = 3, area = 1E5, resid_depth = 0.2,
                 empty_int = 100)

tsd <- abm(730, 1, man_pars = man_pars,
           add_pars = )

plot(VS_conc ~ time, data = tsd, type = 'l')
plot(VSS_conc ~ time, data = tsd, type = 'l')
plot(isCOD_conc ~ time, data = tsd, type = 'l')
plot(sCOD_conc ~ time, data = tsd, type = 'l')
plot(TS_conc ~ time, data = tsd, type = 'l')
tail(tsd)

plot(slurry_mass ~ time, data = tsd, type = 'o')
plot(slurry_depth ~ time, data = tsd, type = 'l')
plot(loss_cum_f_COD ~ time, data = tsd, type = 'l')
plot(loss_cum_f_COD ~ time, data = tsd, type = 'l')
plot(CH4_emis_rate ~ time, data = tsd, type = 'l')

x <- subset(tsd, is.nan(CH4_emis_rate))
tsd$time

matplot(tsd$time, tsd[, nn <- c('m1', 'm2', 'm3')],
  type = 'l', lty = 1, xlab = 'Time (d)', ylab = 'Microbial biomass (g)')
legend('topleft', nn, col = 1:6, lty = 1)

names(tsd)

plot(slurry_mass ~ time, data = tsd2, type = 'l')
plot(COD_conc ~ time, data = tsd, type = 'l')
abline(v = 1:10 * 35, lty = 2)
plot(dsCOD_conc ~ time, data = tsd, type = 'l')
plot(dpCOD_conc ~ time, data = tsd, type = 'l')
plot(ipCOD_conc ~ time, data = tsd, type = 'l')
plot(isCOD_conc ~ time, data = tsd, type = 'l')
plot(m1_conc ~ time, data = tsd, type = 'l')
abline(v = 1:10 * 35, lty = 2)

plot(m2_conc ~ time, data = tsd, type = 'l', xlim = c(0, 730))
lines(tsd2$time + 365, tsd2$m2_conc, type = 'l', col = 'blue')

plot(COD_conc ~ time, data = tsd, type = 'l', xlim = c(0, 730))
lines(tsd2$time + 365, tsd2$COD_conc, type = 'l', col = 'blue')

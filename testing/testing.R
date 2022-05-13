
rm(list = ls())
ff <- list.files('../R', full.names = TRUE)
for (i in ff) source(i)

# Normal abm_regular
out <- abm(100, 1)

out <- abm(100, 1, add_pars = list(slurry_prod_rate = 1000, slurry_rem_rate = 990, slurry_mass = 1000))

head(out, 2)

plot(slurry_mass ~ time, data = out, type = 'l')
plot(slurry_mass ~ time, data = out2, type = 'l')
plot(COD_conc ~ time, data = out, type = 'l')
abline(v = 1:10 * 35, lty = 2)
plot(dsCOD_conc ~ time, data = out, type = 'l')
plot(dpCOD_conc ~ time, data = out, type = 'l')
plot(ipCOD_conc ~ time, data = out, type = 'l')
plot(isCOD_conc ~ time, data = out, type = 'l')
plot(m1_conc ~ time, data = out, type = 'l')
abline(v = 1:10 * 35, lty = 2)

plot(m2_conc ~ time, data = out, type = 'l', xlim = c(0, 730))
lines(out2$time + 365, out2$m2_conc, type = 'l', col = 'blue')

plot(COD_conc ~ time, data = out, type = 'l', xlim = c(0, 730))
lines(out2$time + 365, out2$COD_conc, type = 'l', col = 'blue')

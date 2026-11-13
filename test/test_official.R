rm(list = ls())
library(ggplot2)
library(ABM)


sim_official <- abm(add_pars = list(pH = 6.5, conc_fresh.sulfate = 0, slurry_mass = 0))

saveRDS(sim_official, "sim_official.rds")

sim_official <- readRDS("sim_official.rds")
sim_local <- readRDS("sim_local.rds")


all.equal(sim_official, sim_local)

plot(sim_local$time, sim_local$cum_inhib2, ylim = c(0,1))
lines(sim_official$time, sim_official$cum_inhib2)


lines(sim_official$time, sim_official$sulfide_conc)
plot(sim_local$time, sim_local$sulfide_conc, ylim = c(0,0.17))

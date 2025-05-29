# Add concetrations to abm() output

addConcs <- function(dat, pars) {
  
    dat[, paste0(c(pars$grps, pars$subs, 'VFA'), '_conc')] <- dat[, c(pars$grps, pars$subs, 'VFA')] / dat$slurry_mass
    dat[, paste0(c(pars$grps, pars$subs, 'VFA'), '_eff', '_conc')] <- dat[, paste0(c(pars$grps, pars$subs, 'VFA'), '_eff')] / dat$slurry_mass_eff

    return(dat)

}

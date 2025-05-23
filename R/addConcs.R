# Add concetrations to abm() output

addConcs <- function(dat, pars) {
  
    dat[, paste0(c(pars$grps, 'VSd', 'VFA'), '_conc')] <- dat[, c(pars$grps, 'VSd', 'VFA')] / dat$slurry_mass
    dat[, paste0(c(pars$grps, 'VSd', 'VFA'), '_eff', '_conc')] <- dat[, paste0(c(pars$grps, 'VSd', 'VFA'), '_eff')] / dat$slurry_mass

    return(dat)

}

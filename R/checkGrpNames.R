
checkGrpNames <- function(pars) {

  if (!all.equal(names(pars$yield), names(pars$xa_fresh), names(pars$xa_init), names(pars$decay_rate),
                 names(pars$ks_coefficient), names(pars$resid_enrich), names(pars$qhat_opt), 
                 names(pars$T_opt), names(pars$T_min), names(pars$T_max), 
                 names(pars$ki_NH3_min), names(pars$ki_NH3_max), 
                 names(pars$ki_NH4_min), names(pars$ki_NH4_max), 
                 names(pars$pH_lwr), names(pars$upr))) {
    stop('Microbial groups, i.e., element names in `grp_pars`, must match.')
  }

 }

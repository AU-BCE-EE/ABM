# Names of the (non-derivative) output vector returned by rates_cpp()
# Order must match the order in rates_cpp.cpp
# Microbial group variables are named by group, e.g., rut.m1 (becomes rut_m1 in abm() output)

outNames <- function(pars) {

  grps <- names(pars$qhat_opt)
  v <- function(nm) paste0(nm, '.', grps)

  c('COD_load_rate', 'C_load_rate', 'N_load_rate', 'CH4_emis_rate', 'CO2_emis_rate',
    'H2S_emis_rate', 'NH3_emis_rate_pit', 'NH3_emis_rate_floor', 'N2O_emis_rate',
    paste0('conc_fresh.', colnames(pars$conc_fresh_tab$y)),
    v('xa_fresh'),
    'slurry_prod_rate', 'rain', 'evap', 'sum_rut',
    v('rut'), v('xa'), v('qhat'),
    v('qhat_opt'),
    v('cum_inhib'), v('HAC_inhib'), v('NH3_inhib'), v('NH4_inhib'), v('H2S_inhib'), v('sat_index'),
    'VFA_supply_VFA_demand')

}

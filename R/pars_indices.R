# Order of pars elements that are used in rates_cpp()
# By doing so we can access parms in rates_cpp by indices, which is faster than by names
# Note: rates_cpp() uses positions in this vector, so only add new elements at the end!
rates_pars_names <- c('slurry_prod_rate','t_run','n_mic','temp_standard', 'A', 'E', 'scale', 'R', 'alpha_opt_scale_type', 'alpha_opt_scale_CP',
  'T_opt','T_min','T_max','qhat_opt','decay_rate','g_NH4','pH_inhib_overrule','pH_floor','pH_LL','pH_UL','ki_NH3_min',
  'ki_NH3_max','ki_NH4_min','ki_NH4_max','ki_HAC','ki_H2S_slope','ki_H2S_int','ki_H2S_min','IC50_low','area','floor_area',
  'resp','kl','i_meth','i_sr', 'ks_coefficient','ks_SO4','km_urea','conc_fresh', 'carb','pro','lip','ace', 'hyd', 'ace_sr', 'hyd_sr',
  'carb_resp', 'pro_resp', 'lip_resp','COD_conv','xa_fresh','yield', 'rain', 'evap', 'N2O_emis_rate', 'OM','OM_resp','t_min_graze','t_max_graze','slurry_graze_reduc',
  'temp_C_tab', 'pH_tab', 'conc_fresh_tab', 'xa_fresh_tab', 'out_names', 'mol_COD')

pars_indices <- function(pars){

  p_idx <- match(rates_pars_names, names(pars))

  if (anyNA(p_idx)) {
    stop('Missing parameters: ', paste(rates_pars_names[is.na(p_idx)], collapse = ', '))
  }
  
  # subtract 1 to convert to C++ indexing. 
  p_idx <- p_idx - 1

  return(p_idx)
  
}

# Prepare parameter list for rates_cpp()
# Adds elements that must exist before indices are determined, makes sure numeric vectors are stored as double
# (rates_cpp() reads them directly), and adds parameter indices
prepRatesPars <- function(pars) {

  pars$t_run <- 0
  pars$out_names <- outNames(pars)

  for (i in rates_pars_names) {
    if (is.numeric(pars[[i]]) && !is.double(pars[[i]])) {
      storage.mode(pars[[i]]) <- 'double'
    }
  }

  pars$p_idx <- pars_indices(pars)

  return(pars)

}

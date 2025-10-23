pars_indices <- function(pars){

  # manually define the order of pars elements that are used in rates_cpp.
  # by doing so we can access parms in rates_cpp by indices, which is faster than by names

  pars_elem_ordered <- c('slurry_prod_rate','t_run','n_mic','temp_standard', 'A', 'E', 'scale', 'R', 'alpha_opt_scale_type', 'alpha_opt_scale_CP',
  'T_opt','T_min','T_max','qhat_opt','decay_rate','g_NH4','pH_inhib_overrule','pH_floor','pH_LL','pH_UL','ki_NH3_min',
  'ki_NH3_max','ki_NH4_min','ki_NH4_max','ki_HAC','ki_H2S_slope','ki_H2S_int','ki_H2S_min','IC50_low','area','floor_area',
  'resp','kl','i_meth','i_sr', 'ks_coefficient','ks_SO4','km_urea','conc_fresh', 'carb','pro','lip','ace', 'hyd', 'ace_sr', 'hyd_sr',
  'carb_resp', 'pro_resp', 'lip_resp','COD_conv','xa_fresh','yield', 'rain', 'evap', 'N2O_emis_rate', 'OM','OM_resp')
  
  p_idx <- vector(mode = "numeric", length = length(pars_elem_ordered)) 
  p_idx <- match(pars_elem_ordered, names(pars))
  
  # subtract 1 to convert to C++ indexing. 
  p_idx <- p_idx - 1

  return(p_idx)
  
}